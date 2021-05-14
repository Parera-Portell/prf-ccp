/*
 * ccp.c
 * 
 * Copyright 2021 Joan Antoni Parera Portell
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sacio.h>

#define MAX 5000
#define PI 3.14159265359
#define DEG2RAD 0.017453292520
#define RAD2DEG 57.295779513
#define DEG2KM 111.195
#define KM2DEG 0.00899321
#define EARTHRAD 6371.0


/* Funció per calcular la distància sobre el cercle màxim entre dos
 * punts, així com l'azimut */
void garc(float lat0, float lon0, float lat1, float lon1, float *dist,
		float *az){
			
	float dlon, angle, az_angle, x, y;
	
	lat0 *= DEG2RAD; lon0 *= DEG2RAD; lat1 *= DEG2RAD; lon1 *= DEG2RAD;
	dlon = lon1-lon0;
	
	/* Distància */
	angle = acos(sin(lat0)*sin(lat1)+cos(lat0)*cos(lat1)*cos(dlon));
	*dist = angle*RAD2DEG*DEG2KM;
	
	/* Azimut */
	x = sin(dlon)*cos(lat1);
	y = cos(lat0)*sin(lat1)-sin(lat0)*cos(lat1)*cos(dlon);
	az_angle = atan2(x,y);
	*az	= az_angle*RAD2DEG;
	
	if (*az < 0){
		*az += 360.0;}
}

/* Funció per calcular el punt final a partir d'un punt inicial, un
 * azimut (graus) i una distància (km) sobre el cercle màxim */
void inter(float lat0, float lon0, float az, float dist, float *lat1, 
			float *lon1){
			
	float angle, x, y;
	
	lat0 *= DEG2RAD; lon0 *= DEG2RAD; az *= DEG2RAD; 
	dist *= KM2DEG*DEG2RAD;
	
	angle = asin(sin(lat0)*cos(dist)+cos(lat0)*sin(dist)*cos(az));
	*lat1 = angle*RAD2DEG;
	x = sin(az)*sin(dist)*cos(lat0);
	y = cos(dist)-sin(lat0)*sin(angle);
	angle = lon0+atan2(x,y);
	*lon1 = angle*RAD2DEG;
}

/* Funció per projectar un punt sobre un perfil, en funció de la distància
 * del punt a projectar al punt inicial del perfil (dist) en km, de 
 * l'azimut del punt a projectar al punt inicial del perfil (az0) i de 
 * l'azimut del perfil (az1), en graus. Retorna el punt x sobre el perfil 
 * i la distància y de l'estació respecte del perfil en km */
void proj(float dist, float az0, float az1, float *x, float *y){
	
	float azdiff;
	
	azdiff = az1-az0;
	if (azdiff < 0){
		azdiff += 360.0;}
	azdiff *= DEG2RAD;
	*x = dist*cos(azdiff);
	*y = dist*sin(azdiff);
}


/* Funció principal */
int main(int argc, char **argv){
	
	int t0, tf, fzp, w, n, wsum, size, x0, nlen, nerr, nlay, ncols, 
	max=MAX, i;
	float delta, dx, dz, dep, depmax, beg, p, array[MAX], z, mvp, mvs,
	inilat, inilon, finlat, finlon, len, azim, stla, stlo, baz, hw;
	char rf[300], outdir[300], model[300], pvar[20];
	char *prm = argv[1];
	char *rflist = argv[2];
	FILE *prm_file, *list_file, *mod_file;
	
	if(argc < 3){
		printf("\nUsage: mig3d [par file] [rf list]\n");
		exit(1);
	}
	
	/* --------------Check if the parameter file exists-------------- */
	prm_file = fopen(prm, "r");
	if(prm_file == NULL){
		printf("\nParameter file doesn't exist.\n");
		exit(1);
	}
	
	printf("\n******************* CCP STACKING ********************\n");
	printf("Par. file: %s\n", prm);
	printf("RF file: %s\n", rflist);	
	
	
	/* -------------Read parameter file and store values------------- */
	printf("\n-Parameters-\n");	
	w=0;
	wsum=0;
	w=fscanf(prm_file, "%f", &inilat); wsum += w;
	w=fscanf(prm_file, "%f", &inilon); wsum += w;
	w=fscanf(prm_file, "%f", &finlat); wsum += w;
	w=fscanf(prm_file, "%f", &finlon); wsum += w;
	w=fscanf(prm_file, "%d,%d", &t0, &tf); wsum += w;
	w=fscanf(prm_file, "%f,%f", &dx, &dz); wsum += w;
	w=fscanf(prm_file, "%f", &depmax); wsum += w;
	w=fscanf(prm_file, "%f", &hw); wsum += w;
	w=fscanf(prm_file, "%d", &fzp); wsum += w;
	w=fscanf(prm_file, "%s", outdir); wsum += w;
	w=fscanf(prm_file, "%s", model); wsum += w;
	w=fscanf(prm_file, "%s", pvar); wsum += w;
	fclose(prm_file);
	
	/* Check if all parameters are read */
	if(wsum != 14){
		printf("Error reading parameter file. Exiting...\n");
		exit(1);}
	
	printf("Initial lat/lon: \t%.2f %.2f deg\n", inilat, inilon);
	printf("Final lat/lon: \t\t%.2f %.2f deg\n", finlat, finlon);
	printf("Time window: \t\t%d %d s\n", t0, tf);
	printf("Delta x and delta z: \t%.2f %.2f km\n", dx, dz);
	printf("Maximum depth: \t\t%.2f km\n", depmax);
	printf("Half width: \t\t%.2f km\n", hw);
	printf("Fresnel zone period: \t%d\n", fzp);
	printf("Output: \t\t%s\n", outdir);
	printf("Earth model: \t\t%s\n", model);
	printf("Ray param. variable: \t%s\n", pvar);
	
	/* Check if list file exists */
	list_file = fopen(rflist, "r");
	if(list_file == NULL){
		printf("\nList file doesn't exist.\n");
		exit(1);}
	fclose(list_file);
	
	/* Check if model file exists */
	mod_file = fopen(model, "r");
	if(mod_file == NULL){
		printf("\nModel file doesn't exist.\n");
		exit(1);}
	
	
	/* ------Read and transform gradient model to layered model------ */
	/* Structures to store the Earth model */
	nlay=depmax/dz; /* number of layers */
	struct learthmodel{
			float vp[nlay], vs[nlay];};
	struct earthmodel{
			float dep0[1000], vp0[1000], vs0[1000];};
	struct learthmodel emod; /* model nou */
	struct earthmodel emod0; /* model original */
	
	mod_file = fopen(model, "r");
	w=3;
	n=0;
	/* The model is read into a structure */
	while(w == 3 && dep < depmax){
		w=fscanf(mod_file, "%f,%f,%f", &emod0.dep0[n],&emod0.vp0[n],&emod0.vs0[n]);
		dep=*(emod0.dep0+n);
		n++;
	}
	fclose(mod_file);
	/* Now the gradient model is transformed to a layered model */
	n=0;
	for(w=0; w<nlay; w++){
		z=dz*w;
		i=0; mvp=0; mvs=0;
		/* This loop checks the model dz and compares it to the selected
		 * dz. Also takes into consideration discontinuities (dz==0). */
		while(*(emod0.dep0+n) < z+dz){
			if(n==0){
				mvp+=*(emod0.vp0+n);
				mvs+=*(emod0.vs0+n);
				i++;
			}
			else{
				if(*(emod0.dep0+n) != *(emod0.dep0+n+1)){
					mvp+=*(emod0.vp0+n);
					mvs+=*(emod0.vs0+n);
					i++;
				}
				else{n++;}
			}
			n++;
		}
		/* Assign average velocities if model dz < dz */
		if(i!=0){
			*(emod.vp+w)=mvp/i;
			*(emod.vs+w)=mvs/i;
		}
		/* Assign velocities if model dz >= dz */
		else{n--;
			 *(emod.vp+w)=*(emod0.vp0+n);
			 *(emod.vs+w)=*(emod0.vs0+n);
		}
	}
		/*printf("%f,%f,%f\n",z,*(emod.vp+w),*(emod.vs+w));*/
	
	
	/* ------------------Inicialització del perfil------------------- */
	garc(inilat,inilon,finlat,finlon,&len,&azim);
	ncols = len/dx;
	printf("\nProfile length: %.2f km\n", len);
	printf("Profile azimuth: %.2f deg\n", azim);
	printf("Profile columns: %d\n", ncols);
	float perfil[ncols][nlay];
	
	/* --------------Start looping through Rfs in list--------------- */	
	list_file = fopen(rflist, "r");
	w=1;
	while(w == 1){
		w=fscanf(list_file, "%s", rf);
		
		/* Call rsac1 (SAC library) to read sac file. Returns the array
		 * variable. nlen: array length; beg: beggining time; del: delta
		 * or time sampling; mx: MAX; nerr: error return flag; strlen(file):
		 * length of file path. */
		rsac1(rf, array, &nlen, &beg, &delta, &max, &nerr, strlen(rf));

		/* Check the error status (0=success) */
		if (nerr != 0){
			printf("\nError reading SAC file: %s\n", rf);
			exit (nerr);}
		
		/* Call getfhv (SAC library) to get variables from header */
		getfhv(pvar, &p, &nerr, strlen(pvar));
		getfhv("BAZ", &baz, &nerr, strlen("BAZ"));
		getfhv("STLA", &stla, &nerr, strlen("STLA"));
		getfhv("STLO", &stlo, &nerr, strlen("STLO"));
	
			
		}
	fclose(list_file);
	
	return 0;
}

