/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* A program for calculating linkage disequilibrium stats */
/* Copyright 2005 Daniel Rios & Tim Cutts                 */

#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

#define SIZE 150
#define WINDOW_SIZE 100000
#define INITIAL_LIST_SIZE 256
#define USER_ERROR 1
#define SYSTEM_ERROR 2

/* Macros for fetching particular haplotypes from the allele_counters */
#define AABB allele_counters[0x0000]
#define AABb allele_counters[0x0001]
#define AAbb allele_counters[0x0003]
#define AaBB allele_counters[0x0004]
#define AaBb allele_counters[0x0005]
#define Aabb allele_counters[0x0007]
#define aaBB allele_counters[0x000c]
#define aaBb allele_counters[0x000d]
#define aabb allele_counters[0x000f]

/* Macro which turns pairs of characters into a two bit integer */
#define genotype2int(g) (((g[0] & 0x20) >> 4) | ((g[1] & 0x20) >> 5))

typedef struct {
  int person_id;
  uint8_t genotype;
} Genotype;

typedef struct{
  int locus;
  int position;
  int variation_feature_id;
  int seq_region_end;
  int population_id;
  int number_genotypes;
  Genotype genotypes[SIZE];
} Locus_info;

typedef struct{
  double D;
  double r2;
  double theta;
  int N;
  double d_prime;
  int people;
} Stats;

typedef struct {
  int head;
  int tail;
  int sz;
  Locus_info *locus;
} Locus_list;

typedef struct{
  int number_haplotypes;
  uint8_t haplotype[SIZE];
} Haplotype;

void init_locus_list(Locus_list *l) {
  l->sz = INITIAL_LIST_SIZE;
  l->tail = -1;
  l->head = 0;
  l->locus = malloc(INITIAL_LIST_SIZE*sizeof(Locus_info));
  if (l->locus == NULL) {
    perror("Could not allocate memory");
    exit(SYSTEM_ERROR);
  }
}

void reallocate_locus_list(Locus_list *l) {
  Locus_info *t;
  l->sz *= 2;
  if (( t = realloc(l->locus, l->sz * sizeof(Locus_info))) == NULL) {
    perror("Out of memory reallocating locus list");
    exit(SYSTEM_ERROR);
  }
  l->locus = t;
}

/* dequeue is so simple, it's just a macro */
#define dequeue(ll) ll.head++

void enqueue(Locus_list *ll, int locus, int pos, int var_ft_id, int seq_reg_end, int pop_id, int personid, uint8_t genotype) {
  Locus_info *l;

  ll->tail++;
  if (ll->tail == ll->sz)
    reallocate_locus_list(ll);

  l = &ll->locus[ll->tail];
  l->locus = locus;
  l->position = pos;
  l->variation_feature_id = var_ft_id;
  l->seq_region_end = seq_reg_end;
  l->population_id = pop_id;
  l->genotypes[0].person_id = personid;
  l->genotypes[0].genotype = genotype;
  l->number_genotypes = 1;
}

void major_freqs(const Haplotype * haplotypes, double *f_A, double *f_B){
  int f_a = 0, f_b = 0;
  int f_A_tmp = 0, f_B_tmp = 0;
  int total = 0;
  int i;
  int tmp;
  uint8_t h;

  for (i=0;i<haplotypes->number_haplotypes;i++){
    h = haplotypes->haplotype[i];
    tmp = ((h & 0x8) >> 3) + ((h & 0x4) >> 2);
    f_a += tmp;
    f_A_tmp += (2 - tmp);
    tmp = ((h & 0x2) >> 1) + (h & 0x1);
    f_b += tmp;
    f_B_tmp += (2 - tmp);
    total = total + 2;
  }
  if (total == 0) {
    *f_A = 0.0;
    *f_B = 0.0;
    return;
  }
  if (f_a > f_A_tmp){
    tmp = f_a;
    f_a = f_A_tmp;
    f_A_tmp = tmp;
  }
  if (f_b > f_B_tmp){
    tmp = f_b;
    f_b = f_B_tmp;
    f_B_tmp = f_b;
  }
  *f_A = ((double)f_A_tmp / (double)total);
  *f_B = ((double)f_B_tmp / (double)total);
  return;  
}

int by_person_id(const void *v1, const void *v2){
  Genotype * data1 = (Genotype *)v1;
  Genotype * data2 = (Genotype *)v2;
  
  if (data1->person_id > data2->person_id) return 1;
  if (data1->person_id == data2->person_id) return 0;
  if (data1->person_id < data2->person_id) return -1;
  return 0;
}

void calculate_pairwise_stats(Locus_info *first, Locus_info *second, Stats *s){

  int allele_counters[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  int nAB = 0, nab = 0, nAb = 0, naB = 0;
  int N;
  double theta = 0.5;
  double thetaprev = 2.0;
  double tmp;
  double f_A, f_B;
  double D,r2, Dmax = 0.0, d_prime;
  int i;
  int j;
  Genotype *genotype;
  uint8_t haplotype;

  Haplotype haplotypes;

  /*sort the arrays with the person_id numbers*/
  qsort(first->genotypes,first->number_genotypes,sizeof(Genotype),by_person_id);
  qsort(second->genotypes,second->number_genotypes,sizeof(Genotype),by_person_id);

  for (i=0;i< first->number_genotypes;i++){
    for (j=i;j<second->number_genotypes;j++){
      if (first->genotypes[i].person_id == second->genotypes[j].person_id){
	genotype = &second->genotypes[j];
	/*the second locus has the same person_id*/
	haplotype = first->genotypes[i].genotype << 2;
	haplotype |= genotype->genotype;
	
	allele_counters[haplotype]++;
	
	haplotypes.haplotype[i] = haplotype;
	haplotypes.number_haplotypes = i+1;
	break;
      }    
      else{
	/*they have different, */
	if (first->genotypes[i].person_id < second->genotypes[j].person_id){
	  break;
	}
	/*      fprintf(stderr, "Could not find %hd in locus %d\n",first->genotypes[i].person_id,second->locus);*/
      }
    }
  }

  nAB = 2*AABB + AaBB + AABb;
  nab = 2*aabb + Aabb + aaBb;
  nAb = 2*AAbb + Aabb + AABb;
  naB = 2*aaBB + AaBB + aaBb;
  
  N = nAB + nab + nAb + naB + 2*AaBb;

  /* Initialise stats variable */
  bzero(s, sizeof(Stats));

  if (N < 40){
    /*not enough individuals, return */
    s->N = N;
    return;
  }
  
  while(fabs(theta-thetaprev) > 0.0001){
    thetaprev = theta;
    tmp = ((nAB + (1-theta)*AaBb)*(nab + (1-theta)*AaBb) + (nAb + theta*AaBb)*(naB + theta*AaBb));
    theta = (tmp==0) ? 0.5 : ((nAb + theta*AaBb)*(naB + theta*AaBb))/ tmp;
  }

  /*now calculate stats*/
  major_freqs(&haplotypes,&f_A,&f_B);
  D = (nAB+(1-theta)*AaBb) / N - (f_A*f_B);

  tmp = (f_A*f_B*(1-f_A)*(1-f_B));

  r2 = (tmp==0) ? 0 : D*D/tmp;

  if (D < 0){
    if (f_A*f_B < ((1-f_A)*(1-f_B))) Dmax = f_A*f_B;
    if (f_A*f_B >= ((1-f_A)*(1-f_B))) Dmax = (1-f_A)*(1-f_B);
  }

  if (D >0){
    if (f_A*(1-f_B) < (1-f_A)*f_B) Dmax =  f_A*(1-f_B);
    if (f_A*(1-f_B) >= (1-f_A)*f_B) Dmax = (1-f_A)*f_B;
  }

  s->D = D;
  s->r2 = r2;
  s->theta = theta;
  s->N = N;
  s->d_prime = (Dmax == 0) ? 0.0 : D/Dmax;
  s->people = haplotypes.number_haplotypes;
}

void calculate_ld(const Locus_list *ll, int seq_region_id, FILE *fh){
  Locus_info *next, *head;
  Stats stats;

  /* Doesn't look like it, but sets head and next to the first and
     second entries - I love C, sometimes */

  next = &ll->locus[ll->head];
  head = next++;

  int i;
  for (i = ll->head; i < ll->tail; i++, next++) {
    calculate_pairwise_stats(head, next, &stats);
    if (stats.r2 < 0.05 || stats.N < 40 || stats.r2 > 1)
      continue;
    fprintf(fh, "%d\t%d\t%hd\t%d\t%d\t%d\t%f\t%f\t%d\n",
	    head->variation_feature_id,
	    next->variation_feature_id,
	    head->population_id,
	    seq_region_id,
	    head->position,
	    next->position,
	    stats.r2,
	    fabs(stats.d_prime),
	    stats.N);   
  }
}

int main(int argc, char *argv[])
{
  FILE *fh;
  FILE *in;

  int seq_region_id_previous = -1;
  int locus;
  int position;
  int personid;
  int variation_feature_id;
  int seq_region_id;
  int seq_region_end;
  int population_id;
  char genotype[2];

  Locus_list locus_list;
  Locus_info* l_tmp;

  int last_position = -1;
  int file_line = 1;
  int fields;

  /* Check the command line parameters */
  if ((argc != 3) && (argc != 1)) {
    fprintf(stderr, "Usage:  %s [input_file output_file]\n", argv[0]);
    return USER_ERROR;
  }
  if (argc > 1) {
    if ((in = fopen(argv[1],"r"))==NULL) {
      perror("Could not open input file");
      return USER_ERROR;
    }
    if ((fh = fopen(argv[2],"w"))==NULL) {
      perror("Could not open output file");
      return USER_ERROR;
    }
  } else {
    in = stdin;
    fh = stdout;
  }

  init_locus_list(&locus_list);
  
 loop: 
  while (!feof( in )){
    if (8 != (fields = fscanf(in,
                              "%d%d%d\t%2c\t%d%d%d%d",
                              &locus,
                              &position,
                              &personid,
                              genotype,
                              &variation_feature_id,
                              &seq_region_id,
                              &seq_region_end,
                              &population_id))) {
        if (fields == EOF)
            break;

        fprintf(stderr, "Input parse failure at line %d got %d fields\n",
                file_line, fields);
        return USER_ERROR;
    }
    file_line++;
    
    /* Check both are 'a' or 'A' */
    if ((genotype[0] | genotype[1] | 0x20) != 'a') {
      fprintf(stderr, "Genotype must be AA, Aa or aa, not %s\n",genotype);
      return USER_ERROR;
    }

    /* Make all hets the same order */
    if (genotype[0] == 'a' && genotype[1] == 'A') {
      genotype[0] = 'A'; genotype[1] = 'a';
    }
    
    /*new region, calculate the ld information for the remaining markers*/
    if (seq_region_id != seq_region_id_previous){
      while (locus_list.tail >= locus_list.head){
	calculate_ld(&locus_list, seq_region_id_previous, fh);
	dequeue(locus_list);
      }
      seq_region_id_previous = seq_region_id;
      last_position = -1;
      /* Can reset the queue to the beginning */
      locus_list.head = 0;
      locus_list.tail = -1;
    }

    /*still the same marker, add the genotypes to the array and get next one*/
    if (position == last_position){
      l_tmp = &locus_list.locus[locus_list.tail];
      l_tmp->genotypes[l_tmp->number_genotypes].person_id = personid;
      l_tmp->genotypes[l_tmp->number_genotypes].genotype = genotype2int(genotype);
      l_tmp->number_genotypes++;
      if (l_tmp->number_genotypes == SIZE) {
          fprintf(stderr, "Number of genotypes supported by the program (%d) exceeded\n", SIZE);
          return SYSTEM_ERROR;
      }
      goto loop;
    }

    /*check if the new position is farther than the limit.*/
    /*if so, calculate the ld information for the values in the array*/
    while ((locus_list.tail >= locus_list.head) &&
	    (abs(locus_list.locus[locus_list.head].position - position) > WINDOW_SIZE)){
      calculate_ld(&locus_list, seq_region_id_previous, fh);
      dequeue(locus_list);  
    }
    if (locus_list.tail < locus_list.head) {
      /* Can reset the queue to the beginning */
      locus_list.head = 0;
      locus_list.tail = -1;
    }

    /*new marker*/
    last_position = position;
    
    enqueue(&locus_list, locus, position, variation_feature_id,
	    seq_region_end, population_id, personid, genotype2int(genotype));

  }

  if (in != stdin)
    fclose(in);

  while(locus_list.tail >= locus_list.head){
    calculate_ld(&locus_list, seq_region_id_previous, fh);
    dequeue(locus_list);
  }

  if (fh != stdout)
    fclose(fh);

  return 0;
}
