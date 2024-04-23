/***** Single population *****/
#include<math.h>
#include<iostream>
#include<string.h>
#include<stdio.h>
#include<stdlib.h> 
#include<iomanip> 				
#include<fstream>
#include<tuple>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <unistd.h>

using namespace std;

#define PI 3.1415927

#define NIL (0)    

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long idum;
gsl_rng *gerador;


struct populacao{
  int N;
  int N0;
  double mut;
  int *step;
  int nstep;
  int *contind;
  int ntraits;
  double sigma_r;
  double sigma_trait;
  double sigma_trait_mut;
  int n_subpopul;
  double alpha_r;
  int K;
  double migration;
  double lambda;
};

struct caminhada{

  int nstep;
  double fitness_final;
  int sequence_final;
};

struct otimo{
  
  double *traits;
  int tau;
  double *v;
  double *trait_shift;
  double phenotype;
  
};



struct sequencia{
  //  int *cont_mut;
  int *ind;
  int L;
  int ind_min;
  int *global_phase;
  int ind_max;
  int ind_max_ant;
  int ind_max_landscape;
  int dham;
  int Nmax;
  int *max;
  double roughness;
  double rough_local;
  //  int *path;
  int path_size;
  double *fitness;
  double **traits;
  double **mut_trait;
  double W_max;
  double *trait_init;
  int *number;
  int **offspring;
  int **parent;
  int *z;
  int **C;
  double fmin;
  char **seq;
  int **epistasis;
  double **f;
  double *phenotype;
  double range;
  int K;
  int cont_mut;
};

void fitnessinit(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum);
unsigned int pot(int a1, int a2);
void phenotype(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum);
void dynamics(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum);
double Err(double aa,int bb);
double fitness(int ind, struct sequencia *sequence);
double ran2();


int i1=0;

int main(int ac, char **av)
{
  FILE *ptt1, *ptt2, *ptt3;

  otimo optimum;

  sequencia sequence;

  populacao popul;

  caminhada walk;

  int nwalks, max1, cont_path_total, cont_path, cont_newpath, i, j, k, cont_walk, verif, *nstep_path, *predictability, sum, **step_path, sum_step, max, conf, config, t, tmax, model_dominant, *cont_seq, nstep_real, int_conf, min_step, cont_max, *step_real, n_landscapes, pos, oper1, ind1, *count_ending, *pilha_ending, n_ending, indice_ending, popul_max, parallelism, parallelism_ending, indice, ind_sub, ind_max, *s_path, **path_descent, end_point, end_point_1, ind_path, j_sup, accessibility, cont_extinction, cont_nonextinction, maxsize, cont_g, conf_land;

  double fit_max, fit_med, soma1, P2, *path_divergence, med_step, *W, sum_P2, sum_P2_2, cont_conf, sum_step_land, sum_step_land2, med_P2, med_P2_2, var_P2, erro_P2, var_step, error_step, diverg, med_diverg, med_diverg_2, var_diverg, error_diverg, diverg_sym, med_diverg_sym, med_diverg_sym_2, var_diverg_sym, error_diverg_sym, sum_path_divergence_sym, sum_path_divergence_sym_2, sum_path_divergence, sum_path_divergence_2, sum_trait_2, drop_fitness, *trait, sum_t_2, sum_divergence, sum_parallelism, sum_parallelism_ending, sum_land_parallelism, sum_land_parallelism_ending, sum_land_divergence, sum_land_divergence_2, sum_land_P2, sum_land_P2_2, med_parallelism, med_parallelism_ending, med_divergence, med_divergence_2, var_divergence, error_divergence, error_P2, distance_phen, divergence, f_max, med_accessibility, sum_accessibility, sum_land_P2_path, sum_land_P2_2_path, P2_path, sum_land_step, sum_land_step_2, med_P2_path, med_P2_2_path, var_P2_path, error_P2_path,med_step_2, sum_cont_extinction, sum_divergence_2, sum_step_2, inc_phenotype, r;

  long seed;  //define a seed

  char arq1[1000], arq2[1000], arq3[1000];

  if (ac!=11)
    {
      cout  <<  "start the program like this:\n" << av[0]
            << " <N> <U> <s> <sd> <alpha> <Nlevel> <tmax> <Tterm> <config> <Nbins> <Nmax> <semente> \n"
            << endl;
      exit (-1);
    }

  j=0;
  
  sequence.L = atoi (av[++j]);
  sequence.K =  atoi (av[++j]);
  popul.K = atoi (av[++j]);
  popul.mut = atof (av[++j]);
  n_landscapes = atoi (av[++j]);
  conf_land = atoi (av[++j]);
  popul.alpha_r = atof (av[++j]);
  drop_fitness = atof (av[++j]);
  sequence.W_max = atof (av[++j]);
  sequence.range = atof (av[++j]);
    
  cout  << "#invocation: ";
  for (int i=0; i<ac; i++){
    cout << av[i] << " ";
  }
  cout << endl;

  seed = time (NULL) * getpid();  //set the seed to system time
  gerador=gsl_rng_alloc(gsl_rng_mt19937);  
  gsl_rng_set(gerador, seed); //give the seed to random generator
  
  srand(time(NULL));

  max1 = (int) pow(2.,(double)(sequence.L));
  if( sequence.L==32 )
    max1 = RAND_MAX;

  max = (int)pow(2.,(double)(sequence.K+1));

  sum_trait_2 = -2*log((1-drop_fitness)/1.);

  sequence.epistasis =  new int*[sequence.L];
  for( i=0; i<sequence.L; i++ )
      sequence.epistasis[i] =  new int[sequence.K+1];

  sequence.f =  new double*[sequence.L];
  for( i=0; i<sequence.L; i++ )
    sequence.f[i] =  new double[max];

  sequence.number = new int[max1];

  sequence.phenotype = new double[max1];

  sequence.fitness = new double[max1];

  sequence.ind = new int[100000];

  popul.step =  new int[10000];

  step_real =  new int[10000];

  nstep_path = new int[1000];

  step_path = new int*[1000];
  for( i=0; i<1000; i++ )
    step_path[i] =  new int[500];

  predictability = new int[10000];

  pilha_ending = new int[10000];

  count_ending = new int[10000];
  
  W = new double[10000];

  config = n_landscapes;
  conf = 0;
  cont_conf = 0;

  sum_cont_extinction = 0;

  cont_extinction = 0;
  cont_nonextinction = 0;
  n_ending = 0;
  accessibility = 0;
  sum_step = 0;
  while( (++conf)<=config )
    {
      popul.N = popul.K;    

      if( (conf%conf_land)==1 )
	phenotype(&sequence,&popul,&optimum);

      sequence.ind_min = (int)(gsl_ran_flat(gerador, 0., 1.)*max1);

      inc_phenotype =  pow((-2*log((1-drop_fitness)/sequence.W_max)),0.5);

      optimum.phenotype = sequence.phenotype[sequence.ind_min] + inc_phenotype;
  
      f_max = 0;
      for( indice=0; indice<max1; indice++ )
	{
	  distance_phen = (sequence.phenotype[indice]-optimum.phenotype)*(sequence.phenotype[indice]-optimum.phenotype);
	  
	  sequence.fitness[indice] = sequence.W_max*exp(-distance_phen/(2*popul.alpha_r*popul.alpha_r));
	  if( sequence.fitness[indice]>f_max )
	    {
	      f_max = sequence.fitness[indice];
	      ind_max = indice;
	    }
	  
	}

      fitnessinit(&sequence,&popul,&optimum);
      
      t = 0;
      model_dominant = sequence.ind_min;
      popul.nstep = 0;
      popul.step[0] = 0;
      while( (popul.N!=0) )
	{

	  dynamics(&sequence,&popul,&optimum);
	  
	  if( (popul.N==0) )
	    {
	      cont_extinction++;
	      break;
	    }

	  cont_g = 0;
	  for( i=0; i<popul.N; i++ )
	    if( sequence.fitness[sequence.ind[i]]>1 )
	      cont_g++;

	  if( cont_g>200 )
	    {
	      cont_nonextinction++;
	      break;
	    }
	  
	}

      if( ind_max==sequence.ind_min )
	conf--;

    }

 
  sprintf(arq2,"Rescue-NK-L%d-K%d-N%d-U%g-Wmax%g-range%g.dat",sequence.L,sequence.K,popul.K,popul.mut,sequence.W_max,sequence.range);
  ptt2 = fopen(arq2,"a");
  fprintf(ptt2,"%g \t  %g  \n",drop_fitness,((double)cont_extinction/n_landscapes));
  fclose(ptt2);

}


void phenotype(sequencia *sequence, populacao *popul, otimo *optimum)
{

  int *pilha, i, j, j1, max, aux, k, indice, segmento1, config, max1, oper1, pos;
  double aleat;

  max = (int)pow(2.,(double)(sequence->K+1));
  max1 = (int) pow(2.,(double)sequence->L);

  pilha = new int[sequence->L+1];

  for( i=0; i<sequence->L; i++ )
    {
      for( j=0; j<sequence->L; j++ )
	{
          if( j<i )
	    pilha[j] = j;
          else if( j>i )
	    pilha[j-1] = j;
	}

           
      for( j=0; j<sequence->K; j++ )
	{
	  aux = (int)(ran2()*(sequence->L-j-1));
	  sequence->epistasis[i][j] = pilha[aux];
          for( j1=aux; j1<(sequence->L-j-2); j1++ )
	    pilha[j1] = pilha[j1+1];
	}
      sequence->epistasis[i][sequence->K] = i;  

    }


  for( i=0; i<sequence->L; i++ )
    for( j=0; j<max; j++ )
      sequence->f[i][j] = gsl_ran_flat(gerador, 0., sequence->range);


  for( indice=0; indice<max1; indice++ )
    {
      sequence->phenotype[indice] = 0.;
      for( i=0; i<sequence->L; i++ )
	{
	  config = 0;
	  for( j=0; j<=sequence->K; j++ )
	    {
	      segmento1 = sequence->epistasis[i][j]%(sequence->L);
	      config += ((indice >> segmento1) & 1)*pot(2,j);  
	    }
	  sequence->phenotype[indice] += sequence->f[i][config];
	}
      sequence->phenotype[indice] /= sequence->L;
    }


  delete[] pilha;
  
}




void fitnessinit(sequencia *sequence, populacao *popul, otimo *optimum)
{
  int i, max, max1, ind, j;

  max1 = (int) pow(2.,(double)(sequence->L));
  if( sequence->L==32 )
    max1 = RAND_MAX;

  popul->N = popul->K;
  for( j=0; j<popul->N; j++ )
    sequence->ind[j] = sequence->ind_min; 

  
}


void dynamics(sequencia *sequence, populacao *popul, otimo *optimum)
{
  int i, k, m, contbin, kinf, aleat, dig, oper, j, bit, verif, k1, cont, nlabelaux, size, indtau, *ndel, auxind, *ind1, pos, oper1, ind_pilha, max1, nmut, n_offspring, *pilha, seq_ind;

  unsigned int *N_newborn;
     
  double x, soma, aux, soma1, aux1, mut, Fit, *auxF, malp, mutb, mb, r, Fmax, fitness, saux, F_tot, distance_2, med, a, factor;

  double *prob;

  max1 = (int) pow(2.,(double)(sequence->L));
  if( sequence->L==32 )
    max1 = RAND_MAX;

  ind1 = new int[1000000];
  
  pilha = new int[100000];

  for( i=0; i<max1; i++ )
    sequence->number[i] = 0;

  for( i=0; i<popul->N; i++ )
    sequence->number[sequence->ind[i]]++;

  ind_pilha = 0;
  for( i=0; i<max1; i++ )
    if( sequence->number[i]>0 )
      {
	pilha[ind_pilha] = i;
	ind_pilha++;
      }
  
  cont = 0;
  for( i=0; i<ind_pilha; i++ )
    {
      seq_ind = pilha[i];
      fitness = sequence->fitness[seq_ind]*sequence->number[seq_ind];
      //     cout << fitness << "\t" << seq_ind << "\t" << sequence->fitness[seq_ind] << endl;
      n_offspring = gsl_ran_poisson(gerador,fitness);
      for( j=0; j<n_offspring; j++ )
	{
	  ind1[cont] = seq_ind;
	  r = gsl_ran_flat(gerador, 0., 1.);
	  if( r < popul->mut )
	    {
	      pos = (int)(gsl_ran_flat(gerador, 0., 1.)*sequence->L);
	      oper1 = pot(2,pos);
	      ind1[cont] = (ind1[cont]^oper1);
	      //  sequence->cont_mut++;
	    }
	  cont++;
	}
    }

  //  cout << cont << endl;

  popul->N = cont;
  for( i=0; i<popul->N; i++ )
    sequence->ind[i] = ind1[i];

  delete[] ind1;

  delete[] pilha;
 
}





double fitness(int ind, sequencia *sequence)
{
  int k, config, indice, segmento1;
  double fitness;

  fitness = sequence->fitness[ind];

  return(fitness);
}



double Err(double aa,int bb)
{
  double erro;
  
  erro = pow((double)aa/bb,0.5);
  
  return(erro);
}




unsigned int hamming_distance(unsigned int sequence1, unsigned  int sequence2)
{
	unsigned int dist = 0;
	unsigned int val  = sequence1^sequence2;
	
	while(val)
	{
		++dist;
		val &= (val-1);
	}
	
	return dist;
}


double path_divergence(int* path1, int npath1, int* path2, int npath2) 
{
    int distance = 0;
    cout << npath1 << "\t" << npath2 << endl;
    for(int i=0; i<npath1; i++)
    {
        int min_dist = hamming_distance(path1[i], path2[0]);
        for(int j=0; j<npath2; j++)
        {
            if(min_dist == 0) break;
            
            int dist = hamming_distance(path1[i], path2[j]);
            min_dist = (min_dist < dist)?min_dist:dist;
        }
        distance += min_dist;
    }
    return double(distance)/npath1;
}

double path_divergence_sym(int* path1, int npath1, int* path2, int npath2) 
{
    int distance;
    int distance1 = 0;
    int distance2 = 0;
    
    for(int i=0; i<npath1; i++)
    {
        int min_dist1 = hamming_distance(path1[i], path2[0]);
        for(int j=0; j<npath2; j++)
        {
            if(min_dist1 == 0) break;
            
            int dist1 = hamming_distance(path1[i], path2[j]);
            min_dist1 = (min_dist1 < dist1)?min_dist1:dist1;
        }
        distance1 += min_dist1;
    }
    
    for(int i=0; i<npath2; i++)
    {
        int min_dist2 = hamming_distance(path2[i], path1[0]);
        for(int j=0; j<npath1; j++)
        {
            if(min_dist2 == 0) break;
            
            int dist2 = hamming_distance(path2[i], path1[j]);
            min_dist2 = (min_dist2 < dist2)?min_dist2:dist2;
        }
        distance2 += min_dist2;
    }

    distance = distance1 + distance2;
    return double(distance)/(npath1+npath2);
}




unsigned int pot(int a1, int a2)
{
   int i;

   unsigned int produto;

   produto = 1;
   for( i=1; i<=a2; i++ )
      produto = produto*a1;
   return(produto);
}

double ran2()
{

         int j;
         long k;
         static long idum2=123456789;
         static long iy=0;
         static long iv[NTAB];
         float temp;

	 //	 cout << idum << endl; 
         if (idum <= 0) {
                 if (-(idum) < 1) idum=1;
                 else idum = -(idum);
                 idum2=(idum);
                 for (j=NTAB+7; j>=0; j--) {
                         k=(idum)/IQ1;
                         idum=IA1*(idum-k*IQ1)-k*IR1;
                         if (idum < 0) idum += IM1;
                         if (j < NTAB) iv[j] = idum;

                 }
                 iy=iv[0];// cout << idum<< endl;
         }
         k=(idum)/IQ1;
         idum=IA1*(idum-k*IQ1)-k*IR1;
         if (idum < 0) idum += IM1;
         k=idum2/IQ2;
         idum2=IA2*(idum2-k*IQ2)-k*IR2;
         if (idum2 < 0) idum2 += IM2;
         j=iy/NDIV;
         iy=iv[j]-idum2;
         iv[j]=idum;
         if (iy < 1) iy += IMM1;
         if ((temp=AM*iy) > RNMX) return RNMX;
         else return temp;
}

