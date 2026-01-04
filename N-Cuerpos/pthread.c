#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

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

cuerpo_t *cuerpos;
float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
float **fuerzaX_por_hilo, **fuerzaY_por_hilo, **fuerzaZ_por_hilo;
int N, pasos, P;
int dt = 1.0f;
float toroide_alfa;
float toroide_theta;
float toroide_incremento;
float toroide_lado;
float toroide_r;
float toroide_R;

pthread_barrier_t barrera;

double dwalltime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

void calcularFuerzas_paralelo(int idHilo) {
    int i, j;
    float dif_X, dif_Y, dif_Z;
    float distancia, F;

//cada hilo trabaja con un cuerpo distinto (hilo 0: cuerpo 0, cuerpo 4, cuerpo 8, etc.)
    for (i = idHilo; i < N - 1; i += P) {
	//calcula la fuerza de ese cuerpo con todos los cuerpos siguientes.
        for (j = i + 1; j < N; j++) {
		//si estan en la misma posicion se omite el calculo.
            if (cuerpos[i].px == cuerpos[j].px &&
                cuerpos[i].py == cuerpos[j].py &&
                cuerpos[i].pz == cuerpos[j].pz) continue;
		
		//calculamos la distancia entre cuerpos
            dif_X = cuerpos[j].px - cuerpos[i].px;
            dif_Y = cuerpos[j].py - cuerpos[i].py;
            dif_Z = cuerpos[j].pz - cuerpos[i].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);
							
		//calculamos la fuerza resultante
            F = (G * cuerpos[i].masa * cuerpos[j].masa) / (distancia * distancia);
		
            dif_X *= F;
            dif_Y *= F;
			dif_Z *= F;

		//para cada hilo guardamos la fuerza calculada de cada cuerpo.
            fuerzaX_por_hilo[idHilo][i] += dif_X;
            fuerzaY_por_hilo[idHilo][i] += dif_Y;
			fuerzaZ_por_hilo[idHilo][i] += dif_Z;

            fuerzaX_por_hilo[idHilo][j] -= dif_X;
            fuerzaY_por_hilo[idHilo][j] -= dif_Y;
			fuerzaZ_por_hilo[idHilo][i] -= dif_Z;
        }
    }
}

void moverCuerpos_paralelo(int idHilos) {
    float ax, ay;
		

	//cada hilo trabaja con un cuerpo distinto (hilo 0: cuerpo 0, cuerpo 4, cuerpo 8, etc.) 
    for (int c = idHilos; c < N; c += P) {
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
        ax = fuerza_totalX[c] / cuerpos[c].masa;
        ay = fuerza_totalY[c] / cuerpos[c].masa;
		//az = fuerza_totalZ[c] / cuerpos[c].masa;
	
	//calculamos la velocidad
        cuerpos[c].vx += ax * dt;
        cuerpos[c].vy += ay * dt;
		//cuerpos[c].vz += az * dt;
	
	//calculamos la posicion
        cuerpos[c].px += cuerpos[c].vx * dt;
        cuerpos[c].py += cuerpos[c].vy * dt;
		//cuerpos[c].pz += cuerpos[c].vz * dt;

        fuerza_totalX[c] = 0.0;
        fuerza_totalY[c] = 0.0;
		fuerza_totalZ[c] = 0.0;
    }
}

void* worker(void* arg) {
    int idHilo = *(int*)arg;

    for (int paso = 0; paso < pasos; paso++) {      
	
	//calculamos fuerzas
        calcularFuerzas_paralelo(idHilo);
	
	//esperamos que todos los hilos terminen
        pthread_barrier_wait(&barrera);

	//calculamos el movimiento de los cuerpos.
        moverCuerpos_paralelo(idHilo);

	//esperamos a que todos los hilos terminen.
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

        fuerza_totalX[cuerpo] = 0.0;
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;

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

	cuerpos = (cuerpo_t*)malloc(sizeof(cuerpo_t) * N);
    fuerza_totalX = (float*)malloc(sizeof(float) * N);
    fuerza_totalY = (float*)malloc(sizeof(float) * N);
    fuerza_totalZ = (float*)malloc(sizeof(float) * N);

//matriz de fuerzas
// Reservar memoria para las filas, cantidad de workers
     fuerzaX_por_hilo= (float**)malloc(sizeof(float*) * P);
     fuerzaY_por_hilo= (float**)malloc(sizeof(float*) * P);
     fuerzaZ_por_hilo= (float**)malloc(sizeof(float*) * P);
    // Reservar memoria para las columnas, cada cuerpo
    for(int i = 0; i < P; i++) {
    fuerzaX_por_hilo[i] = (float*)malloc(N* sizeof(float));
	fuerzaY_por_hilo[i] = (float*)malloc(N * sizeof(float));
	fuerzaZ_por_hilo[i] = (float*)malloc(N * sizeof(float));
    }

	for (int i = 0; i < N; i++) {
		for(int idHilo=0; idHilo<P; idHilo++){
            fuerzaX_por_hilo[idHilo][i] = 0.0;
            fuerzaY_por_hilo[idHilo][i] = 0.0;
	    	fuerzaZ_por_hilo[idHilo][i] = 0.0;
        }
	}
}

void finalizar(){
	
	 free(cuerpos);
    free(fuerza_totalX);
    free(fuerza_totalY);
    free(fuerza_totalZ);
    for(int i = 0; i < P; i++){
	free(fuerzaX_por_hilo[i]);
	free(fuerzaY_por_hilo[i]);
	free(fuerzaZ_por_hilo[i]);	
	}
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

	inicializar();

    inicializarCuerpos(cuerpos,N);

    pthread_barrier_init(&barrera, NULL, P);

    pthread_t hilos[P];
    int ids[P];

    double tInicio = dwalltime();

    for (int i = 0; i < P; i++) {
        ids[i] = i;
        pthread_create(&hilos[i], NULL, &worker,(void*)&ids[i]);
    }

    for (int i = 0; i < P; i++) {
        pthread_join(hilos[i], NULL);
    }

    double tFin = dwalltime();
    printf("Tiempo total: %f segundos\n", tFin - tInicio);

	
    pthread_barrier_destroy(&barrera);
   
	finalizar();
    return 0;
}




























