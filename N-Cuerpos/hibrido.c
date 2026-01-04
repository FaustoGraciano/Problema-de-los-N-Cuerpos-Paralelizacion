#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <mpi.h>
#include <stddef.h>
#include <string.h>

#define G 6.673e-11
#define PI (3.141592653589793)
#define ESTRELLA 0
#define POLVO 1
#define H2 2 //Hidrogeno molecular

typedef struct cuerpo {
    float masa;
    float px, py, pz;
    float vx, vy, vz;
    float r, g, b;
    int cuerpo;
} cuerpo_t;

cuerpo_t *cuerpos_locales;
cuerpo_t *cuerpos_remoto;

float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
float **fuerzaX_por_hilo, **fuerzaY_por_hilo, **fuerzaZ_por_hilo;

int N, pasos, P;
float dt = 1.0f;

float toroide_alfa;
float toroide_theta;
float toroide_incremento;
float toroide_lado;
float toroide_r;
float toroide_R;

///HIBRIDO
typedef struct {
	float x,y,z;
}tf;

int blockSize;
int rank;
int size;

tf **fuerza_por_hilo_remota;
tf *fuerzas_acumuladas;


pthread_barrier_t barrera;

double dwalltime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

//calculo de fuerzas entre mi bloque y un bloque remoto
void calcularFuerzas_paralelo_remoto(int idHilo, cuerpo_t *cuerpos_remoto) {
    int i, j;
    float dif_X, dif_Y, dif_Z;
    float distancia, F;

    //cada hilo trabaja con un cuerpo distinto (hilo 0: cuerpo 0, cuerpo 4, cuerpo 8, etc.)
    for (i = idHilo; i < blockSize; i += P) {
        //calcula la fuerza de ese cuerpo con todos los cuerpos siguientes.
        for (j = 0; j < blockSize; j++) {
            //si estan en la misma posicion se omite el calculo.
            if (cuerpos_locales[i].px == cuerpos_remoto[j].px &&
                cuerpos_locales[i].py == cuerpos_remoto[j].py &&
                cuerpos_locales[i].pz == cuerpos_remoto[j].pz) continue;

            //calculamos la distancia entre cuerpos
            dif_X = cuerpos_remoto[j].px - cuerpos_locales[i].px;
            dif_Y = cuerpos_remoto[j].py - cuerpos_locales[i].py;
            dif_Z = cuerpos_remoto[j].pz - cuerpos_locales[i].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

            //calculamos la fuerza resultante
            F = (G * cuerpos_locales[i].masa * cuerpos_remoto[j].masa) / (distancia * distancia);

            dif_X *= F;
            dif_Y *= F;
			dif_Z *= F;

            //para cada hilo guardamos la fuerza calculada de cada cuerpo.
            fuerzaX_por_hilo[idHilo][i] += dif_X;
            fuerzaY_por_hilo[idHilo][i] += dif_Y;
			fuerzaZ_por_hilo[idHilo][i]	+= dif_Z;
			
			fuerza_por_hilo_remota[idHilo][j].x -= dif_X;
			fuerza_por_hilo_remota[idHilo][j].y -= dif_Y;
			fuerza_por_hilo_remota[idHilo][j].z -= dif_Z;
			

        }
    }
}


/////////////////////////////


void calcularFuerzas_paralelo(int idHilo) {
    int i, j;
    float dif_X, dif_Y, dif_Z;
    float distancia, F;

//cada hilo trabaja con un cuerpo distinto (hilo 0: cuerpo 0, cuerpo 4, cuerpo 8, etc.)
    for (i = idHilo; i < blockSize - 1; i += P) { 
	//calcula la fuerza de ese cuerpo con todos los cuerpos siguientes.
        for (j = i + 1; j < blockSize; j++) {
		//si estan en la misma posicion se omite el calculo.
            if (cuerpos_locales[i].px == cuerpos_locales[j].px &&
                cuerpos_locales[i].py == cuerpos_locales[j].py &&
                cuerpos_locales[i].pz == cuerpos_locales[j].pz) continue;
		
		//calculamos la distancia entre cuerpos
            dif_X = cuerpos_locales[j].px - cuerpos_locales[i].px;
            dif_Y = cuerpos_locales[j].py - cuerpos_locales[i].py;
            dif_Z = cuerpos_locales[j].pz - cuerpos_locales[i].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);
							
		//calculamos la fuerza resultante
            F = (G * cuerpos_locales[i].masa * cuerpos_locales[j].masa) / (distancia * distancia);
		
            dif_X *= F;
            dif_Y *= F;
			dif_Z *= F;

		//para cada hilo guardamos la fuerza calculada de cada cuerpo.
            fuerzaX_por_hilo[idHilo][i] += dif_X;
            fuerzaY_por_hilo[idHilo][i] += dif_Y;
			fuerzaZ_por_hilo[idHilo][i] += dif_Z;		

            fuerzaX_por_hilo[idHilo][j] -= dif_X;
            fuerzaY_por_hilo[idHilo][j] -= dif_Y;
			fuerzaZ_por_hilo[idHilo][j] -= dif_Z;
        }
    }
}

void moverCuerpos_paralelo(int idHilos) {
    float ax, ay;
	
	//cada hilo trabaja con un cuerpo distinto (hilo 0: cuerpo 0, cuerpo 4, cuerpo 8, etc.) 
    for (int c = idHilos; c < blockSize; c += P) {
		//por cada cuerpo se suma la fuerza de cada hilo
		for (int hilo = 0; hilo < P; hilo++) {
        	fuerza_totalX[c] += fuerzaX_por_hilo[hilo][c];
            fuerza_totalY[c] += fuerzaY_por_hilo[hilo][c];
		    fuerza_totalZ[c] += fuerzaZ_por_hilo[hilo][c];
			
			fuerzaX_por_hilo[hilo][c] = 0.0;
            fuerzaY_por_hilo[hilo][c] = 0.0;
	    	fuerzaZ_por_hilo[hilo][c] = 0.0;
        }
       	
	//calculamos la aceleracion 
        ax = fuerza_totalX[c] / cuerpos_locales[c].masa;
        ay = fuerza_totalY[c] / cuerpos_locales[c].masa;
		//az = fuerza_totalZ[c] / cuerpos_locales[c].masa;
	
	//calculamos la velocidad
        cuerpos_locales[c].vx += ax * dt;
        cuerpos_locales[c].vy += ay * dt;
		//cuerpos[c].vz += az * dt;
	
	//calculamos la posicion
        cuerpos_locales[c].px += cuerpos_locales[c].vx * dt;
        cuerpos_locales[c].py += cuerpos_locales[c].vy * dt;
		//cuerpos_locales[c].pz += cuerpos_locales[c].vz * dt;

        fuerza_totalX[c] = 0.0;
        fuerza_totalY[c] = 0.0;
		fuerza_totalZ[c] = 0.0;
    }
}

void* worker(void* arg) {
    int idHilo = *(int*)arg;

    for (int paso = 0; paso < pasos; paso++) {
       
		 // Paso 1: enviar mis cuerpos a todos los procesos con menor rank
        if (idHilo == 0) {
            for (int otro_rank = 0; otro_rank < rank; otro_rank++) {
                MPI_Send(cuerpos_locales, blockSize * sizeof(cuerpo_t), MPI_BYTE, otro_rank, 99, MPI_COMM_WORLD);
            }
        }
		
        //calcular fuerzas internas
        calcularFuerzas_paralelo(idHilo);
        pthread_barrier_wait(&barrera);

        

        // Paso 2: recibir cuerpos remoto de procesos con mayor rank, calcular fuerzas remotas, enviar fuerzas calculadas
        for (int otro_rank = rank + 1; otro_rank < size; otro_rank++) {
            if (idHilo == 0) {
                MPI_Recv(cuerpos_remoto, blockSize * sizeof(cuerpo_t), MPI_BYTE, otro_rank, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            pthread_barrier_wait(&barrera);
			
			//calculo fuerzas con bloques remotos.
            calcularFuerzas_paralelo_remoto(idHilo, cuerpos_remoto);
            pthread_barrier_wait(&barrera);

            if (idHilo == 0) {		//creo fuerzas acumuladas para no pasar matriz.
			for (int c = 0; c < blockSize; c++) {
                    fuerzas_acumuladas[c].x = 0.0f;
                    fuerzas_acumuladas[c].y = 0.0f;
					fuerzas_acumuladas[c].z = 0.0f;
                    for (int h = 0; h < P; h++) {
                        fuerzas_acumuladas[c].x += fuerza_por_hilo_remota[h][c].x;
                        fuerzas_acumuladas[c].y += fuerza_por_hilo_remota[h][c].y;
						fuerzas_acumuladas[c].z += fuerza_por_hilo_remota[h][c].z;

                        // Limpiar para la próxima iteración
                        fuerza_por_hilo_remota[h][c].x = 0.0f;
                        fuerza_por_hilo_remota[h][c].y = 0.0f;
						fuerza_por_hilo_remota[h][c].z = 0.0f;
                    }
                }
                MPI_Send(fuerzas_acumuladas, 2 * blockSize, MPI_FLOAT, otro_rank, 100, MPI_COMM_WORLD);
			}

        }

        // Paso 3: recibir fuerzas remotas de procesos con menor rank
        for (int otro_rank = 0; otro_rank < rank; otro_rank++) {
             if(idHilo == 0) {
                MPI_Recv(fuerzas_acumuladas, 2 * blockSize, MPI_FLOAT, otro_rank, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            pthread_barrier_wait(&barrera);

            for (int c = idHilo; c < blockSize; c += P) {
                fuerza_totalX[c] += fuerzas_acumuladas[c].x;
                fuerza_totalY[c] += fuerzas_acumuladas[c].y;
				fuerza_totalZ[c] += fuerzas_acumuladas[c].z;
            }
        }

        pthread_barrier_wait(&barrera);

        // Paso 4: mover cuerpos
        moverCuerpos_paralelo(idHilo);
        pthread_barrier_wait(&barrera);
    }

    return NULL;
}




void inicializarEstrella(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa =0.001*8;

        if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
        }else{
            toroide_alfa+=toroide_incremento;
        }

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);

    	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

		cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
		cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
		cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarPolvo(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001*4;
	
        if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
        }else{
            toroide_alfa+=toroide_incremento;
        }

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);
	
	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;
    
	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 0.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 0.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarH2(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001;

	if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
	}else{
            toroide_alfa+=toroide_incremento;
	}

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);

	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarCuerpos(cuerpo_t *cuerpos,int N){
 int cuerpo;
 double n = N;

	toroide_alfa = 0.0;
	toroide_theta = 0.0;
	toroide_lado = sqrt(N);
	toroide_incremento = 2*M_PI / toroide_lado;
	toroide_r = 1.0;
	toroide_R = 2*toroide_r;
	
	srand(time(NULL));

	for(cuerpo = 0; cuerpo < N; cuerpo++){

		cuerpos[cuerpo].cuerpo = (rand() %3);

		if (cuerpos[cuerpo].cuerpo == ESTRELLA){
			inicializarEstrella(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == POLVO){
			inicializarPolvo(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == H2){
			inicializarH2(&cuerpos[cuerpo],cuerpo,n);
		}

	}

		cuerpos[0].masa = 2.0e2;
	        cuerpos[0].px = 0.0;
		cuerpos[0].py = 0.0;
		cuerpos[0].pz = 0.0;
		cuerpos[0].vx = -0.000001;
		cuerpos[0].vy = -0.000001;
		cuerpos[0].vz = 0.0;

		cuerpos[1].masa = 1.0e1;
	        cuerpos[1].px = -1.0;
		cuerpos[1].py = 0.0;
		cuerpos[1].pz = 0.0;
		cuerpos[1].vx = 0.0;
		cuerpos[1].vy = 0.0001;
		cuerpos[1].vz = 0.0;
}

void inicializar(){
	
	cuerpos_remoto = (cuerpo_t*)malloc(sizeof(cuerpo_t) * blockSize);

    cuerpos_locales = (cuerpo_t*)malloc(sizeof(cuerpo_t) * blockSize);
	
	fuerzas_acumuladas = (tf*)malloc(sizeof(tf) * blockSize);

	
    fuerza_totalX = (float*)malloc(sizeof(float) * blockSize);
    fuerza_totalY = (float*)malloc(sizeof(float) * blockSize);
    fuerza_totalZ = (float*)malloc(sizeof(float) * blockSize);

	//matriz de fuerzas
	// Reservar memoria para las filas, cantidad de workers
     fuerzaX_por_hilo= (float**)malloc(sizeof(float*) * P);
     fuerzaY_por_hilo= (float**)malloc(sizeof(float*) * P);
     fuerzaZ_por_hilo= (float**)malloc(sizeof(float*) * P);

	//FUERZAS POR HILO REMOTA
	fuerza_por_hilo_remota= (tf**)malloc(sizeof(tf*) * P);

    // Reservar memoria para las columnas, cada cuerpo
    for(int i = 0; i < P; i++) {
		fuerzaX_por_hilo[i] = (float*)malloc(blockSize * sizeof(float));
		fuerzaY_por_hilo[i] = (float*)malloc(blockSize * sizeof(float));
		fuerzaZ_por_hilo[i] = (float*)malloc(blockSize * sizeof(float));

		fuerza_por_hilo_remota[i] = (tf*)malloc(blockSize * sizeof(tf));
    }

	for (int i = 0; i < blockSize; i++) {
		fuerza_totalX[i] = 0.0;
		fuerza_totalY[i] = 0.0;
		fuerza_totalZ[i] = 0.0;		
		
		for(int idHilo=0; idHilo<P; idHilo++){
            fuerzaX_por_hilo[idHilo][i] = 0.0;
            fuerzaY_por_hilo[idHilo][i] = 0.0;
	    	fuerzaZ_por_hilo[idHilo][i] = 0.0;
			fuerza_por_hilo_remota[idHilo][i].x = 0.0;
			fuerza_por_hilo_remota[idHilo][i].y = 0.0;
			fuerza_por_hilo_remota[idHilo][i].z = 0.0;
        }
	}

}

void finalizar(){
		
	free(cuerpos_remoto);
    free(cuerpos_locales);
	free(fuerzas_acumuladas);
    free(fuerza_totalX);
    free(fuerza_totalY);
    free(fuerza_totalZ);
    for(int i = 0; i < P; i++){
	free(fuerzaX_por_hilo[i]);
	free(fuerzaY_por_hilo[i]);
	free(fuerzaZ_por_hilo[i]);
	free(fuerza_por_hilo_remota[i]);	
	}
	free(fuerza_por_hilo_remota);
    free(fuerzaX_por_hilo);
    free(fuerzaY_por_hilo); 
    free(fuerzaZ_por_hilo);	
}


int main(int argc, char* argv[]) {
    if (argc < 5) {
        printf("Uso: %s <N cuerpos> <DT> <pasos> <hilos>\n", argv[0]);
        return -1;
    }

    N = atoi(argv[1]);
    dt = atoi(argv[2]);
    pasos = atoi(argv[3]);
    P = atoi(argv[4]);
	
	//inicio MPI
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
		
	///hibrido divido cuerpos por size
	blockSize=(N/size);

	
	inicializar();
	
	//si soy rank=0 inicializo todos los cuerpos.
	cuerpo_t *cuerpos_totales = NULL;
	if (rank == 0) {
		cuerpos_totales = malloc(sizeof(cuerpo_t) * N);
		inicializarCuerpos(cuerpos_totales, N);
	}


    pthread_barrier_init(&barrera, NULL, P);

    pthread_t hilos[P];
    int ids[P];
	//arranca el tiempo.
	double tInicio;
	if (rank==0) tInicio = dwalltime();

	//se distribuye los cuerpos entre los ranks.
	MPI_Scatter(cuerpos_totales, blockSize * sizeof(cuerpo_t), MPI_BYTE, cuerpos_locales, blockSize * sizeof(cuerpo_t),MPI_BYTE, 0, MPI_COMM_WORLD);

	//creo hilos
    for (int i = 0; i < P; i++) {
        ids[i] = i;
        pthread_create(&hilos[i], NULL, &worker,(void*)&ids[i]);
    }

    for (int i = 0; i < P; i++) {
        pthread_join(hilos[i], NULL);
    }
	
	//se recibe los bloques calculados de los ranks y los junta en cuerpos_totales.
	MPI_Gather(cuerpos_locales,blockSize * sizeof(cuerpo_t), MPI_BYTE,cuerpos_totales,blockSize * sizeof(cuerpo_t), MPI_BYTE, 0,MPI_COMM_WORLD);


   
	if (rank == 0) {    
		double tFin = dwalltime(); 
		printf("Tiempo total: %f segundos hibrido\n", tFin - tInicio);}

    pthread_barrier_destroy(&barrera);
	
	if (rank == 0) {
    free(cuerpos_totales);
	}

	finalizar();

	MPI_Finalize();
    return 0;
}
