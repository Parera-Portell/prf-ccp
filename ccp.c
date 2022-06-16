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
 * Variables de la funció main:
 * t0: temps inicial; tf: temps final; w-wsum: variables per controlar
 * el nombre de paràmetres que es llegeixen del fitxer de paràmetres;
 * i-j-k-n: variables comptador; nlen: mida de la matriu RF; nerr: variable
 * per control d'errors; nlay-ncols: nombre de files i columnes del
 * perfil; delta: freqüència de les dades en segons; dx-dz: increment de
 * distància lateral i vertical del perfil; depmin-depmax: profunditat
 * mínima i màxima del perfil; dep: variable de control de la profunditat; 
 * beg: temps entre inici de les dades i arribada d'ona P; p: paràmetre 
 * de raig sísmic; z: variable de control de la profunditat del perfil;
 * mvp-mvs: velocitats sísmiques del model; inilat-inilon-finlat-finlon:
 * latitud i longitud d'inici i final del perfil; len: longitud de perfil;
 * azim: azimut de perfil; stla-stla: coordenades de l'estació; baz:
 * backazimuth; hw: mostreig lateral del perfil en km; az0-az1: azimuts
 * per càlculs geogràfics; dist0-dist1: distància per càlculs geogràfics;
 * x-y-xi-yi: punts projectats respecte el perfil en km; lati-loni:
 * coordenades per càlculs geogràfics; ds: distància epicentral entre dos
 * punts al llarg del raig sísmic; dt: temps de viatge entre dos punts al
 * llarg del raig sísmic; a-b-c-l-q: variables auxiliars per càlculs;
 * up-us: velocitats P i S transformades a un model de capes planes;
 * zf: profunditat transformada a un model de capes planes; psph: paràmetre
 * de raig sísmic en s/rad; acctime: temps de viatge acumulat; amp:
 * amplitud de l'RF; fzr: radi de la zona de fresnel; wl: longitud d'ona;
 * iangle: angle d'incidència; xamp: amplitud corregida.
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <sacio.h>
#include <fftw3.h>

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

/* Funció per projectar un punt sobre un perfil, en funció de la distància
 * del punt a projectar al punt inicial del perfil (dist) en km, de 
 * l'azimut del punt inicial del perfil al punt a projectar (az0) i de 
 * l'azimut del perfil (az1), en graus. Retorna el punt x sobre el perfil 
 * i la distància y de l'estació respecte del perfil en km */
void proj(float dist, float az0, float az1, float *x, float *y){
	
	float azdiff;
	
	azdiff = az1-az0;
	azdiff *= DEG2RAD;
	*x = dist*cos(azdiff);
	*y = dist*sin(azdiff);
}


/* Funció principal */
void main(int argc, char **argv){
	
	int t0, w, n, wsum, nlen, nerr, nlay, ncols, max=MAX, i, j, k, fz;
	float delta, dx, dz, dep, depmin, depmax, beg, p, array[MAX], z, mvp, 
	mvs, inilat, inilon, finlat, finlon, len, azim, stla, stlo, baz, hw, 
	az0, az1, dist0, dist1, x, y, xi, yi, ds, dt, a, b, c, up, us, zf, nu, 
	acctime, amp, latj, lonj, fzr, wl, xamp, dzf, stx, l, q;
	char rf[500], outfile[500], outres[500], model[500], pvar[20];
	char *prm = argv[1];
	char *rflist = argv[2];
	FILE *prm_file, *list_file, *mod_file, *out_file, *out_res;
	/* Arrays de l'FFT */
	fftw_complex fftin[MAX], fftout[MAX], phase;
	fftw_plan plan1, plan2;
	
	if(argc < 3){
		printf("\nUsage: mig3d [par file] [rf list]\n");
		exit(1);
	}
	
	/* ---------Comprovar si existeix el fitxer de paràmetres-------- */
	prm_file = fopen(prm, "r");
	if(prm_file == NULL){
		printf("\nParameter file doesn't exist.\n");
		exit(1);
	}
	
	printf("\n*******************************************************\n");
	printf("******************** CCP STACKING *********************\n");
	printf("*******************************************************\n");
	printf("by Joan A. Parera Portell (2021)\n");
	printf("\nPar. file: %s\n", prm);
	printf("RF file: %s\n", rflist);	
	
	
	/* ---------------Lectura del fitxer de paràmetres--------------- */
	printf("\n-Parameters-\n");	
	w=0;
	wsum=0;
	w=fscanf(prm_file, "%f", &inilat); wsum += w;
	w=fscanf(prm_file, "%f", &inilon); wsum += w;
	w=fscanf(prm_file, "%f", &finlat); wsum += w;
	w=fscanf(prm_file, "%f", &finlon); wsum += w;
	w=fscanf(prm_file, "%d", &t0); wsum += w;
	w=fscanf(prm_file, "%f,%f", &dx, &dz); wsum += w;
	w=fscanf(prm_file, "%f,%f", &depmin, &depmax); wsum += w;
	w=fscanf(prm_file, "%f", &hw); wsum += w;
	w=fscanf(prm_file, "%s", outfile); wsum += w;
	w=fscanf(prm_file, "%s", outres); wsum += w;
	w=fscanf(prm_file, "%s", model); wsum += w;
	w=fscanf(prm_file, "%s", pvar); wsum += w;
	w=fscanf(prm_file, "%f", &nu); wsum += w;
	fclose(prm_file);
	
	/* Comprovar si tots els paràmetres s'han llegit correctament */
	if(wsum != 15){
		printf("Error reading parameter file. Exiting...\n");
		exit(1);}
	
	printf("Initial lat/lon: \t%.2f %.2f deg\n", inilat, inilon);
	printf("Final lat/lon: \t\t%.2f %.2f deg\n", finlat, finlon);
	printf("Initial time: \t\t%d s\n", t0);
	printf("Delta x and delta z: \t%.2f %.2f km\n", dx, dz);
	printf("Min/max depths: \t%.2f %.2f km\n", depmin, depmax);
	printf("Half width: \t\t%.2f km\n", hw);
	printf("Output slice: \t\t%s\n", outfile);
	printf("Sampling slice: \t%s\n", outres);
	printf("Earth model: \t\t%s\n", model);
	printf("Ray param. variable: \t%s\n", pvar);
	printf("Phase weight exp. term: %.2f\n", nu);
	
	/* Comprovar si la llista de RFs existeix */
	list_file = fopen(rflist, "r");
	if(list_file == NULL){
		printf("\nList file doesn't exist.\n");
		exit(1);}
	fclose(list_file);
	
	/* Comprovar si el fitxer del model existeix */
	mod_file = fopen(model, "r");
	if(mod_file == NULL){
		printf("\nModel file doesn't exist.\n");
		exit(1);}
	
	/* Lectura i transformació del model de gradient a model de capes */
	nlay=depmax/dz; /* nombre de capes */
	/* Definició d'estructures per desar el model */
	struct learthmodel{
			float vp[nlay], vs[nlay];};
	struct earthmodel{
			float dep0[1000], vp0[1000], vs0[1000];};
	struct learthmodel emod; /* model nou */
	struct earthmodel emod0; /* model original */
	
	mod_file = fopen(model, "r");
	w=3;
	n=0;
	/* Es desa el model a la memòria */
	while(w==3){
		w=fscanf(mod_file, "%f,%f,%f", &emod0.dep0[n],&emod0.vp0[n],&emod0.vs0[n]);
		n++;
	}
	fclose(mod_file);
	/* Comença la transformació a un model de capes */
	n=0;
	for(w=0; w<nlay; w++){
		z=dz*w;
		i=0; mvp=0; mvs=0;
		
		/* Cas 1: dz al model < dz al perfil */
		while(emod0.dep0[n+1] < z+dz){
			if(emod0.dep0[n] != emod0.dep0[n+1]){
				mvp+=emod0.vp0[n];
				mvs+=emod0.vs0[n];
				i++;
			}
			else{n++;}
			n++;
		}
		/* Assignació de velocitats mitjanes si Cas 1 */
		if(i!=0){
			emod.vp[w]=mvp/i;
			emod.vs[w]=mvs/i;
		}
		/* Cas 2: dz al model >= dz al perfil */
		else{emod.vp[w]=emod0.vp0[n];
			 emod.vs[w]=emod0.vs0[n];
		}
	}

	/* ------------------Inicialització del perfil------------------- */
	garc(inilat,inilon,finlat,finlon,&len,&azim);
	ncols = len/dx;
	printf("\nProfile length: %.2f km\n", len);
	printf("Profile azimuth: %.2f deg\n", azim);
	printf("Profile size: %dx%d\n", nlay, ncols);
	float perfil[ncols][nlay], pnorm[ncols][nlay], absphase[ncols][nlay];
	fftw_complex pphase[ncols][nlay];
	for(i=0; i<nlay; i++){
		for(w=0; w<ncols; w++){
			perfil[w][i] = 0;
			pnorm[w][i] = 1;
		}
	}
	
	/* --------------Start looping through Rfs in list--------------- */
	list_file = fopen(rflist, "r");
	w = 0;
	printf("\nProcessing RFs...");
	while(fgets(rf,500,list_file) != NULL){
		sscanf(rf,"%s",rf);
		/* Call rsac1 (SAC library) to read sac file. Returns the array
		 * variable. nlen: array length; beg: beggining time; del: delta
		 * or time sampling; mx: MAX; nerr: error return flag; strlen(file):
		 * length of file path. */
		rsac1(rf, array, &nlen, &beg, &delta, &max, &nerr, strlen(rf));

		/* Check the error status (0=success) */
		if (nerr != 0){
			printf("\nError reading SAC file: %s\n", rf);
			exit (nerr);}
			
		if(w==0){
			/* FFT directa*/
			plan1 = fftw_plan_dft_1d(nlen, fftin, fftin, FFTW_FORWARD, FFTW_ESTIMATE);
			/* FFT inversa*/
			plan2 = fftw_plan_dft_1d(nlen, fftin, fftout, FFTW_BACKWARD, FFTW_ESTIMATE);}
		w++;
		
		/* Call getfhv (SAC library) to get variables from header */
		getfhv(pvar, &p, &nerr, strlen(pvar));
		getfhv("BAZ", &baz, &nerr, strlen("BAZ"));
		getfhv("STLA", &stla, &nerr, strlen("STLA"));
		getfhv("STLO", &stlo, &nerr, strlen("STLO"));
		
		/* Inicialització de la matriu que sobre la qual es farà l'FFT */
		for(j=0;j<nlen;j++){
			fftin[j] = array[j];
		}
		/* Execució de la FFT directa */
		fftw_execute(plan1);
		/* Càlcul de la senyal analítica: s'eliminen les freqüències
		 * negatives i es deixen a 0. S'assumeix que nlen%2=0 */
		for(j=0;j<nlen;j++){
			if(j>=nlen/2){
				fftin[j] = 0;
			}
		}
		fftw_execute(plan2);
		/* Renormalització de les dades */
		/* Càlcul de la fase instantània (senyal analítica / amplitud
		 * instantània) */
		for(j=0;j<nlen;j++){
			fftout[j] /= (nlen/2);
			fftout[j] /= cabs(fftout[j]);
		}

		/* Punt inicial: coordenades de l'estació projectades al perfil */
		garc(inilat,inilon,stla,stlo,&dist0,&az0);
		proj(dist0,az0,azim,&x,&y);
		stx = x;
		
		/* Inici dels càlculs pel CCCP */
		/* Comprovació si el temps d'inici és positiu o negatiu (respecte
		 * la P directa */
		if(beg<0){acctime = -beg;}
		else{acctime = beg;}
		for(j=0; j<nlay; j++){
			/* Càlcul de la posició del raig sísmic per cada
			 * profunditat del model. A continuació es transformen
			 * els paràmetres a un model de Terra plana */
			/* Increment de profunditat transformada */
			a = -EARTHRAD*log((EARTHRAD-(j+1)*dz)/EARTHRAD);
			zf = -EARTHRAD*log((EARTHRAD-j*dz)/EARTHRAD);
			dzf = a-zf;
			/* Velocitats sísmiques transformades */
			up = (EARTHRAD/(EARTHRAD-j*dz))*emod.vp[j];
			us = (EARTHRAD/(EARTHRAD-j*dz))*emod.vs[j];
			/* Càlcul de la distància epicentral entre els piercing
			 * points del raig sísmic a z i z+1 (Shearer, 2009) */
			a = sqrt(1/(us*us)-p*p);
			ds = p*dzf/a;
			/* Càlcul del temps de viatge entre z i z+1 (Shearer, 2009).
			 * Per trobar el temps de la RF cal restar Ts-Tp, que és el
			 * retard amb que arriba cada conversió Ps */
			b = sqrt(1/(up*up)-p*p);
			/* Comprovació si s'ha arribat al turning point */
			if(a>0 && b>0){
				/* Ps */
				dt = (a-b)*dzf;
				/*  PpPs */
				/*dt = (a+b)*dzf;*/
				/* PpSs */
				/*dt = (a+a)*dzf;*/
				/* Càlcul de l'amplitud dins l'interval de temps */
				acctime += dt;
				a = round(acctime/delta);
				k=a;
				phase = fftout[k];
				amp = array[k];
				/* Càlcul de la projecció sobre el perfil. 
				 * El periode (P) de la FZ == 1. wl = T*v */
				proj(ds,baz,azim,&xi,&yi);
				x += xi;
				y += yi;
				fz = 1;
				wl = fz*us;
				fzr = sqrt(0.5*wl*zf+0.0625*wl*wl);
				/* Factor de saturació amplitud-profunditat */
				q = pow(zf, 1.5);
				/* Se suavitzen les amplituds aplicant una funció
				/* Gaussiana: G(x)=0.5*std*sqrt(2*pi)*-->
				/* -->*exp(-x**2/2*std**2) */
				if(fabs(y)<=hw){
					for(n=0;n<ncols;n++){
						l = n*dx-x;
						a = 0.5*fz*fzr;
						b = -0.5*(l*l)/(a*a);
						if(b < -20){c=0;}
						else{c = exp(b)/a*sqrt(2*3.14159);}
						/* Correcció de l'amplitud */
						xamp = amp*q*c;
						/* Assignació de l'amplitud al perfil i al
						/* perfil de normalització */
						if(xamp != 0){
							perfil[n][j] += xamp;
							pnorm[n][j] += 1;
							pphase[n][j] += phase;
						}
					}				
				}
			}			
		}	
	}
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
	fftw_cleanup();
	printf("%d!\n", w);
	fclose(list_file);
	
	/* Normalització de les amplituds del perfil */
	a = 0;
	for(i=0; i<nlay; i++){
		for(w=0; w<ncols; w++){
			perfil[w][i] /= sqrt(pnorm[w][i]);
			absphase[w][i] = pow(fabs(pphase[w][i])/pnorm[w][i],nu);
			/* Suavitzat per files */
			if(w>0){
				absphase[w][i] = 0.2*absphase[w][i]+(1-0.2)*absphase[w-1][i];
			}
			perfil[w][i] *= absphase[w][i];
			if(fabs(perfil[w][i])>a){
				a = perfil[w][i];
			}
		}
	}
	for(i=0; i<nlay; i++){
		for(w=0; w<ncols; w++){
			perfil[w][i] /= a;
		}
	}
	for(w=0; w<ncols; w++){
		for(i=0; i<nlay; i++){
			/* Suavitzat per columnes */
			if(i>0){
				perfil[w][i] = 0.4*perfil[w][i]+(1-0.4)*perfil[w][i-1];
			}
		}
	}

	/* Escriptura a un fitxer */
	out_file = fopen(outfile, "w");
	out_res = fopen(outres, "w");
	fprintf(out_file, "%s,%s,%s,%s,%s,%s\n", "x","z","a","lat","lon","zdeg");
	fprintf(out_res, "%s,%s,%s\n", "x","z","a");
	for(i=0; i<nlay; i++){
		for(w=0; w<ncols; w++){
			a = w*dx;
			b = i*dz+dz/2;
			stla=inilat+(finlat-inilat)/ncols*w;
			stlo=inilon+(finlon-inilon)/ncols*w;
			c = b*KM2DEG;
			if(i*dz>=depmin){
				fprintf(out_file, "%f,%f,%f,%f,%f,%f\n", a,b,perfil[w][i],stla,stlo,c);
				fprintf(out_res, "%f,%f,%f\n", a,b,pnorm[w][i]);
			}
		}
	}
	fclose(out_file);
	fclose(out_res);
}

