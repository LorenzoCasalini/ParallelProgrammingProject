/** 
  * File:    kernel_p_rotate_molecola_2.cu 
  * 
  * Author:  Lorenzo Casalini
  * Date:     Summer 2020
  * Summary of File: 
  * This file contains a parallel implementation of function Matc Probe Shape. The parallel kernel
  * is the function place_in_best_angle: each thread rotates the fragment in the angle corrisponding to his 
  * thread id and then evaluates the expansion. 
  */ 
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdbool.h> 
#include <assert.h>

typedef struct
{ 
	char name[40];
	int n_atoms;
	int n_bonds;
	double atoms[500];
	int bonds[300];
}molecola;

typedef struct
{    
	int head;
	int tail;
	int elements[500];
}queue;

#define repetitions 10
#define enable_refiniment  false
#define high_precision_step  1
#define low_precision_step  30
#define threshold  0.2 

inline cudaError_t checkCuda(cudaError_t result)
{
	if (result != cudaSuccess) {
		fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
		assert(result == cudaSuccess);
	}
	return result;
}

/**
 * Populates a struct molecola in memory given the name of the text file
 * @param molecola_name 
 * @param m1 empty struct molecola to be inizialized
 */
void create_molecola(char* molecola_name,molecola* m1) {
	FILE *input_file;
	int line_index = 0;
	char* res;
	char line[500];
	int number_of_atoms;
	int number_of_bounds;
	char path[50];
	strcpy(path,"molecules/");
	strcat(path, molecola_name);
	input_file = fopen(path, "r");
	if (input_file == NULL) {
		printf("fopen non funziona\n");
		return;
	}
	res = fgets(line, 100, input_file);
	fgets(line, 100, input_file);
	fgets(line, 100, input_file);
	char* numero = strtok(line, " ");
	number_of_atoms = atoi(numero);
	numero = strtok(NULL, " ");
	number_of_bounds = atoi(numero);
	m1->n_atoms = number_of_atoms;
	m1->n_bonds = number_of_bounds;
	fgets(line, 100, input_file);
	fgets(line, 100, input_file);
	fgets(line, 100, input_file);
	fgets(line, 100, input_file);
	while(1){
		char * token = strtok(line, " ");
		line_index = atoi(token) - 1;
		token = strtok(NULL, " ");
		token = strtok(NULL, " ");
		m1->atoms[3*line_index] = atof(token);
		token = strtok(NULL, " ");
		m1->atoms[3*line_index+1] = atof(token);
		token = strtok(NULL, " ");
		m1->atoms[3*line_index + 2] = atof(token);
		fgets(line,100,input_file);
		if(strncmp(line,"@<TRIPOS>",5)==0){
			break;
		}
		}
	fgets(line, 100, input_file);
    while (strcmp(line, "@<TRIPOS>SUBSTRUCTURE\n") != 0  && res != NULL && strcmp(res,"\n")!=0) {
		char * token = strtok(line, " ");
		line_index = atoi(token) - 1;
		token = strtok(NULL, " ");
		m1->bonds[2*line_index] = atoi(token);
		token = strtok(NULL, " ");
		m1->bonds[2*line_index+1] = atoi(token);
		res = fgets(line, 100, input_file);
		}
	fclose(input_file);
	strcpy(m1->name,molecola_name);
}


/**
 * Checks if a node is present in the queue
 * @param node index of a node
 * @param queue struct queue 
 */
__host__ __device__ bool isPresent(int node, queue* queue){
	for (int i = 0; i < queue->tail; i++) {
		if (queue->elements[i] == node) {
			return true; }
	}
	return false;
}

/**
 * Does a breadth first search on the graph defined by the bonds of the molecola without the one specified by the parameter bond index.
 * It populates queue with the atoms found in the search.
 * The queue must be already inizialized having as first element the first node of the search
 * @param bond_index the bond we must eliminate from the search 
 * @param molecola 
 * @param queue the queue we populate with the adjacent elements, it must be already inizialized   
 */
__host__ __device__ void bfs(int bond_index, molecola*molecola, queue* queue){
	int node = queue->elements[queue->head];
	int n_bonds = molecola -> n_bonds;
	while(queue->head < queue->tail){
		for(int i = 0; i<n_bonds;i++){
		    if(i!=bond_index){
			    if (molecola->bonds[2 * i] == node && !isPresent(molecola->bonds[2 * i + 1], queue)) {
				    queue->elements[queue->tail] = molecola->bonds[2 * i + 1];
				    queue->tail += 1;
			    }
			    else if(molecola->bonds[2 * i+1] == node && !isPresent(molecola->bonds[2 * i], queue)){
				    queue->elements[queue->tail] = molecola->bonds[2 * i];
				    queue->tail += 1;
			    }
		    }
	    }
		queue->head +=1; 
		node = queue -> elements[queue->head];
	}
}

/**
 * It inizialize a queue and then call a bfs to populate the queue with the nodes adjacent to atom. 
 * @param molecola
 * @param queue to be inizialized, bfs will put in this queue the adjacent nodes
 * @param atom 
 * @param bond_index 
 */
__host__ __device__ void find_adjacent_nodes(molecola* molecola, queue* queue, int atom, int bond_index) {
	queue->elements[0] = atom;
	queue->head = 0;
	queue->tail = 1;
	bfs(bond_index, molecola, queue);
}

/** 
 * Given a molecola and a bound, it initializes two queues, in the first queue there is the left atom of the bond, in the second one the right atom. Then it calls a bfs for each queue and if in the queues
 * there are the same elements it is not a rotamer because the two fragment are connected. If in a queue there is only one atom it is not a rotamer.
 * @param bond_index the index of the bond to check
 * @param molecola 
 * @return isRotamer 
 */
bool isRotamer(int bond_index, molecola* molecola) {
	int first_node, second_node;
	bool isRotamer;
	queue q1;
	queue q2;
	first_node = molecola->bonds[2*bond_index];
	second_node = molecola->bonds[2*bond_index+1];
	q1.tail = 1;
	q1.head = 0;
	q1.elements[0] = first_node;
	q2.tail = 1;
	q2.head = 0;
	q2.elements[0] = second_node;
	bfs(bond_index, molecola, &q1);
	bfs(bond_index, molecola, &q2);
	isRotamer = true;
	for (int i = 0; i < q1.tail; i++) {
		for (int j = 0; j < q2.tail; j++) {
			if (q1.elements[i] == q2.elements[j]){
				isRotamer = false;
			}
		}
	}
	if (q1.tail == 1 || q2.tail == 1) {
		isRotamer = false;
	}
	return isRotamer;
}

/**
 * Given a molecola it check every bound if it is a rotamer, then puts the index of the rotamer in a list
 * @param molecola 
 * @param number_of_rotamers a pointer to an int where it puts the number of rotamer it has found 
 * @return rotamer_list the indeces all the rotamers found
 */
int* find_rotamers(molecola* molecola, int* number_of_rotamers) {
	//sempre chiamarlo con un n_rotamers = 0
	int size = molecola->n_bonds;
	bool* x;
	int n_rotamers = 0;
	int* rotamer_list;
	int rotamer_index = 0;
	x = (bool*)malloc(size* sizeof(int));
	for (int i = 0; i < size; i++) {
		if (isRotamer(i, molecola)) { x[i] = true; }
		else { x[i] = false; }
	}
	for (int i = 0; i < size; i++) {
		if (x[i]) {
			n_rotamers += 1;
		}
	}
	rotamer_list = (int*)malloc(n_rotamers * sizeof(int));
	for (int i = 0; i < size; i++) {
		if (x[i]) {
			rotamer_list[rotamer_index] = i;
			rotamer_index += 1;
		}
	}
	free(x);
	*number_of_rotamers = n_rotamers;
	return rotamer_list;
}


__host__ __device__ void normalise(double* x, double* y, double* z) {
	double w = sqrt(*x * *x + *y * *y + *z * *z);
	*x = *x / w;
	*y = *y / w;
	*z = *z / w;
}

/**
 * Given a rotamer it rotates an atom of a given angle around the axis defined by the rotamer
 * @param molecola
 * @param atom to rotate
 * @param rotamer index of the axis we want to rotate around the atom 
 * @param angle 
 */
__host__ __device__ void rotate_atom(molecola*molecola, int atom, int rotamer, int angle) {
	double px, py, pz, p1x, p1y, p1z, p2x, p2y, p2z, rx, ry, rz, qx, qy, qz;
	double tetha = angle*M_PI / 180;
	int rotamer1_index, rotamer2_index;
	double costheta, sintheta;
	px = molecola->atoms[3 * (atom - 1)];
	py = molecola->atoms[3 * (atom - 1) + 1];
	pz = molecola->atoms[3 * (atom - 1) + 2];
	rotamer1_index = molecola->bonds[2 * rotamer];
	rotamer2_index = molecola->bonds[2 * rotamer + 1];
	p1x = molecola->atoms[3 * (rotamer1_index - 1)];
	p1y = molecola->atoms[3 * (rotamer1_index - 1) + 1];
	p1z = molecola->atoms[3 * (rotamer1_index - 1) + 2];
	p2x = molecola->atoms[3 * (rotamer2_index - 1)];
	p2y = molecola->atoms[3 * (rotamer2_index - 1) + 1];
	p2z = molecola->atoms[3 * (rotamer2_index - 1) + 2];
	rx = p2x - p1x;
	ry = p2y - p1y;
	rz = p2z - p1z;
	px = px - p1x;
	py = py - p1y;
	pz = pz - p1z;
	normalise(&rx, &ry, &rz);
	costheta = cos(tetha);
	sintheta = sin(tetha);
	qx = 0;
	qy = 0;
	qz = 0;
	qx += (costheta + (1 - costheta)* rx*rx)*px;
	qx += ((1 - costheta) * rx * ry - rz * sintheta) * py;
	qx += ((1 - costheta) * rx * rz + ry * sintheta) * pz;

	qy += ((1 - costheta) * rx * ry + rz * sintheta) * px;
	qy += (costheta + (1 - costheta) * ry * ry) * py;
	qy += ((1 - costheta) * ry * rz - rx * sintheta) * pz;

	qz += ((1 - costheta) * rx * rz - ry * sintheta) * px;
	qz += ((1 - costheta) * ry * rz + rx * sintheta) * py;
	qz += (costheta + (1 - costheta) * rz * rz) * pz;

	qx += p1x;
	qy += p1y;
	qz += p1z;
	molecola->atoms[3 * (atom - 1)] = qx;
	molecola->atoms[3 * (atom - 1) + 1] = qy;
	molecola->atoms[3 * (atom - 1) + 2] = qz;
}

/**
 * Calculates the distance between two atoms of a molecola
 * @param index_1 first atom 
 * @param index_2 second atom
 * @return distance 
 */
__host__ __device__ double distance(molecola* molecola, int index_1, int index_2) {
	double distance;
	double x1, y1, z1, x2, y2, z2;
	x1 = molecola->atoms[3 * (index_1)];
	y1 = molecola->atoms[3 * (index_1) + 1];
	z1 = molecola->atoms[3 * (index_1) + 2];
	x2 = molecola->atoms[3 * (index_2)];
	y2 = molecola->atoms[3 * (index_2) + 1];
	z2 = molecola->atoms[3 * (index_2) + 2];
	distance = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
	return distance;
}

/**
 * Given a molecola calculates the sum of the distances of all its atoms
 * @param molecola 
 * @return expansion sum of the distances of the atoms in molecola 
 */
__host__ __device__ double measure_expansion(molecola* molecola) {
	double expansion = 0;
	for (int i = 0; i < molecola->n_atoms; i++) {
		for (int j = 0; j < molecola->n_atoms; j++) {
			if (j > i) {
				expansion += distance(molecola, i, j);
			}
		}
	}
	return expansion;
}

/**
 * Given a rotamer and the index of an belonging to that rotamer, first finds the atoms connected the first one and then rotates them
 * @param molecola
 * @param rotamer the index of the bound which we want to rotate the fragmen around
 * @param atom the index of the atom, belonging to the rotamer, that indicates which part of the molecola we must rotate
 * @param angle angle of rotation 
 */
__host__ __device__ void rotate_molecola(molecola* molecola, int rotamer, int atom, int angle) {
	queue q1;
	find_adjacent_nodes(molecola, &q1, atom, rotamer);
	for (int i = 0; i < q1.tail; i++) {
		rotate_atom(molecola, q1.elements[i], rotamer, angle);
	}
}

/**
 * Calculates if a ligand is feasible, a ligand is feasible if all its atoms' distances are > 0.8 A
 * @param molecola
 */
__host__ __device__ bool is_ligand_feasible(molecola* molecola) {
	for (int i = 0; i < molecola->n_atoms; i++) {
		for (int j = 0; j < molecola->n_atoms; j++) {
			if (j > i) {
				if (distance(molecola, i, j) < 0.8) {
					return false;
				}
			}
		}
	}
	return true;
}


/**
 * Calculates the size of the fragment defined by the rotamer index
 * @param molecola
 * @param bond the index of the rotamer
 * @param index is 1 if we evaluate left fragment, 2 right fragment
 * @return size_pct the size of the fragment w.r.t. the size of molecola
 */
double fragment_size(molecola* molecola, int bond, int index) {
	//index 1 if left fragment , 2 right
	queue q1;
	q1.tail = 1;
	q1.head = 0;
	double size_pct;
	if (index == 1) {
		q1.elements[0] = molecola->bonds[2 * bond];
	}
	else if (index == 2) {
		q1.elements[0] = molecola->bonds[2 * bond + 1];
	}
	else {
		printf("Fragment size: Index must be between 1 and 2");
		return 0;
	}
	bfs(bond, molecola, &q1);
	size_pct = (double)q1.tail / molecola->n_atoms;
	return size_pct;
}

/**
 * This kernel creates a molecola in the local memory of each thread and initializes it as the molecola passed
 * as argument. Then initializes to 0 the array expansions passed as argument. At the end each thread rotates 
 * the fragment, defined by bond and atom, of the angle corrisponding to its thread id and writes the expansion
 * in the array expansions
 * @param mol
 * @param bond index of the bond around which we rotate the fragment
 * @param atom index of the atom, is needed to define which part of the molecola we must rotate
 * @param step  step of the rotations
 * @param min_range of rotations
 * @param max_range of rotations
 * @param expansion uninitialized vector to be filled with the expansion corresponding at the angle of 
 * rotation of the given index
 */
__global__ void place_in_best_angle_p(molecola* mol, int bond, int atom, int step, int min_range, int max_range,double* expansions){
	int thread_id = blockIdx.x*blockDim.x + threadIdx.x; 
	molecola mol2;
	mol2 = *mol;
	if(thread_id<360){
		expansions[thread_id] = 0;
	}
	if((thread_id == min_range || thread_id % step == 0) && thread_id < max_range){
		rotate_molecola(&mol2, bond, atom, thread_id);
		if(is_ligand_feasible(&mol2)){
		expansions[thread_id] = measure_expansion(&mol2);
	}
	}
}

/**
 * An helper function which initializes an array of expansion in the unified memory and then calls place_best_angle_p
 * At the end finds the maximum in the vector of expansions and then rotates the molecola in the angle with the 
 * best expansion. 
 * @param mol
 * @param bond index of the bond around which we rotate the fragment
 * @param atom index of the atom, is needed to define which part of the molecola we must rotate
 * @param step  step of the rotations
 * @param min_range of rotations
 * @param max_range of rotations
 */
void place_in_best_angle(molecola* mol, int bond, int atom, int step, int min_range, int max_range){
	double* expansions;
	double best_expansion = 0;
	int best_angle;
	checkCuda(cudaMallocManaged(&expansions, 360*sizeof(double)));
	place_in_best_angle_p<<<12,32>>>(mol, bond, atom, step, min_range, max_range, expansions);
	checkCuda(cudaDeviceSynchronize());
	for(int i = 0; i< 360; i++){
		if(expansions[i]>best_expansion){
			best_expansion = expansions[i];
			best_angle = i;
		}
	}
	checkCuda(cudaFree(expansions));
	rotate_molecola(mol, bond, atom, best_angle);
}

/**
 * Places molecola in the middle of each tile, then evaluate the expansion of the molecola with the fragment placed in the tile, 
 * at the end returns the index of the tile with the best expansion
 * @param mol
 * @param n_tiles number of tile to evaluate
 * @param bond index of the bond around which we rotate the fragment
 * @param atom index of the atom, is needed to define which part of the molecola we must rotate
 * @return best_expansion_tile index of the tile with best expansion
 */
int find_best_tile(molecola* mol, int n_tiles, int bond, int atom) {
	molecola mol2;
	mol2 = *mol;
	int tile_size;
	double expansion;
	double best_expansion=0;
	int best_expansion_tile=0;
	tile_size = floor(360 / n_tiles);
	rotate_molecola(&mol2, bond, atom, tile_size / 2);
	if (is_ligand_feasible(&mol2)) {
		best_expansion = measure_expansion(&mol2);
		best_expansion_tile = 0;
	}
	for (int i = 1; i < n_tiles; i++) {
		rotate_molecola(&mol2, bond, atom, tile_size);
		if (is_ligand_feasible(&mol2)) {
			expansion = measure_expansion(&mol2);
			if (expansion > best_expansion) {
				best_expansion = expansion;
				best_expansion_tile = i;
			}
		}
	}
	return best_expansion_tile;
}

/**
 * Given a molecola, this function find the rotamers and for each rotamer rotates first the left fragment and then the left one, in a way that depends on the fragment size
 * @param molecola a struct molecola
 */
double match_probe_shape(molecola* molecola) {
	int* rotamer_list;
	int rotamer_index = 0;
	int n_rotamers = 0;
	int best_tile;
	int n_tiles = 18;
	int tile_size = 360/n_tiles;
	//non sono sicuro sia giusto
	rotamer_list = find_rotamers(molecola, &n_rotamers);
	for (int j = 0; j < repetitions; j++) {
		for (int i = 0; i < n_rotamers; i++) {
			rotamer_index = rotamer_list[i];
			if (fragment_size(molecola, rotamer_index, 1) < threshold) {
				place_in_best_angle(molecola, rotamer_index, molecola->bonds[2 * rotamer_index], low_precision_step, 0, 360);
			}
			else {
				if (enable_refiniment) {
					best_tile = find_best_tile(molecola, n_tiles, rotamer_index, molecola->bonds[2 * rotamer_index]);
					place_in_best_angle(molecola, rotamer_index, molecola->bonds[2 * rotamer_index], high_precision_step, best_tile*tile_size, (best_tile + 1)*tile_size);
				}
				else {
					place_in_best_angle(molecola, rotamer_index, molecola->bonds[2 * rotamer_index], high_precision_step, 0, 360);
				}
			}
			if (fragment_size(molecola, rotamer_index, 2) < threshold) {
				place_in_best_angle(molecola, rotamer_index, molecola->bonds[2 * rotamer_index + 1], low_precision_step, 0, 360);
			}
			else {
				if (enable_refiniment) {
					best_tile = find_best_tile(molecola, n_tiles, rotamer_index, molecola->bonds[2 * rotamer_index + 1]);
					place_in_best_angle(molecola, rotamer_index, molecola->bonds[2 * rotamer_index + 1], high_precision_step, best_tile*tile_size, (best_tile + 1)*tile_size);
				}
				else {
					place_in_best_angle(molecola, rotamer_index, molecola->bonds[2 * rotamer_index + 1], high_precision_step, 0, 360);
				}
			}
		}
	}
	free(rotamer_list);
	return measure_expansion(molecola);
}
	

int main() {
	clock_t begin = clock();
	int n_molecole = 48;
	double espansion;
	int deviceId;
	cudaGetDevice(&deviceId); 
	molecola* m1;
	//molecola list_of_molecole[1];
	char* molecole_list[] = {"Aspirin.mol2","Diclofenac.mol2","Diplosalsalate.mol2","Flurbiprofen.mol2","Focalin.mol2","Losmiprofen.mol2","Melatonin.mol2","Myfortic.mol2",
	"Nifuradene.mol2","Oxybenzone.mol2","Propiomazine.mol2","Raloxifene.mol2","Relacatib.mol2", "Ribasphere.mol2","Roxoperone.mol2","Sulindac.mol2",
	"1b9v_deposited_1.mol2", "1br6_deposited_1.mol2","1bxq_ligand.mol2", "1c1b_deposited_1.mol2","1ctr_deposited_1.mol2","1cvu_deposited_1.mol2","1cx2_deposited_1.mol2",
    "1ezq_deposited_1.mol2", "1fcx_deposited_1.mol2", "1fl3_deposited_1.mol2", "1fm6_deposited_1.mol2","1fm9_deposited_1.mol2","1fmz_ligand.mol2","1fq5_deposited_1.mol2",
	"1gvx_ligand.mol2", "1gwx_deposited_1.mol2","1h23_ligand.mol2", "1hp0_deposited_1.mol2","1hvy_deposited_1.mol2", "1iiq_ligand.mol2","1lpz_deposited_1.mol2", 
	"1mq6_deposited_1.mol2","1oyt_deposited_1.mol2", "1pso_deposited_1.mol2","1s19_deposited_1.mol2","1uml_deposited_1.mol2","1ydt_deposited_1.mol2","2hnx_ligand.mol2",
	"3l3n_ligand.mol2", "3nhi_ligand.mol2","4djp_ligand.mol2","4gid_ligand.mol2"};
	for (int i = 0; i < n_molecole; i++) {
		checkCuda(cudaMallocManaged(&m1, sizeof(molecola)));
		create_molecola(molecole_list[i],m1);
		espansion = measure_expansion(m1);
		printf("Before expansion Molecola: %s , espansion: %f\n", m1->name,espansion);
		espansion = match_probe_shape(m1);
		printf("Molecola: %s, expansion: %f\n", m1->name, espansion);
		checkCuda(cudaFree(m1));
	}
	clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nTime spent: %f\n", time_spent);
	return 0;
}