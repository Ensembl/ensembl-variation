/* A program for calculating linkage disequilibrium stats */

/*
 Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 Copyright [2016] EMBL-European Bioinformatics Institute
 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at
 
      http://www.apache.org/licenses/LICENSE-2.0
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 How to use
   - Give the command the directories/files you want to search
   - Ignores binary files if possible
   - Code will report its best guess at what type of licence was applied
   - Setting the environment variable APPLY_LICENSE will cause the code to write an Apache2 licence to the file
   - Does not support JavaScript or CSS files
*/


#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>

#include "tbx.h"
#include "vcf.h"

#define SIZE 5000
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
  int position;
  char *var_id;
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

void enqueue(Locus_list *ll, int pos, char *var_id, int pop_id, int personid, uint8_t genotype) {
  Locus_info *l;

  ll->tail++;
  if (ll->tail == ll->sz)
    reallocate_locus_list(ll);

  l = &ll->locus[ll->tail];
  l->position = pos;
  l->var_id = var_id;
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
/*if (f_a > f_A_tmp){
    tmp = f_a;
    f_a = f_A_tmp;
    f_A_tmp = tmp;
  }
  if (f_b > f_B_tmp){
    tmp = f_b;
    f_b = f_B_tmp;
    f_B_tmp = f_b;
  }*/
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
//  fprintf(stderr, "%d\n",  first->number_genotypes);
  int allele_counters[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  int nAB = 0, nab = 0, nAb = 0, naB = 0;
  int N;
  double theta = 0.5;
  double thetaprev = 2.0;
  double tmp;
  double f_A, f_B;
  double D,r2, Dmax = 0.0;
  int i = 0;
  int j = 0;
  int z = 0;
  Genotype *genotype;
  uint8_t haplotype;

  Haplotype haplotypes;

  while (i<first->number_genotypes && j<second->number_genotypes){
    if (first->genotypes[i].person_id == second->genotypes[j].person_id){
      genotype = &second->genotypes[j];
      /*the second locus has the same person_id*/
      haplotype = first->genotypes[i].genotype << 2;
      haplotype |= genotype->genotype;
      
      allele_counters[haplotype]++;
      
      haplotypes.haplotype[z] = haplotype;
      haplotypes.number_haplotypes = z+1;
      z++;
      i++;
      j++;
    }    
    else{
      /*they have different, */
      if (first->genotypes[i].person_id < second->genotypes[j].person_id){
        i++;
      }
      else{
        j++;
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

void calculate_ld(const Locus_list *ll, int seq_region_id, FILE *fh, char *variant){
  Locus_info *next, *head;
  Stats stats;
  
  /* Doesn't look like it, but sets head and next to the first and
     second entries - I love C, sometimes */

  next = &ll->locus[ll->head];
  head = next++;
  int i;
  for (i = ll->head; i < ll->tail; i++, next++) {
    if (variant[0] != 0) {
      if ((strcmp(head->var_id, variant) != 0 ) && (strcmp(next->var_id, variant) != 0))
        continue;
    }
    if (head->population_id != next->population_id)
      continue;
    calculate_pairwise_stats(head, next, &stats);
    if ((float) stats.r2 < 0.05 || stats.N < 40 || (float) stats.r2 > 1 || (float) stats.d_prime > 1)
      continue;
    fprintf(fh, "%d\t%d\t%d\t%s\t%d\t%s\t%f\t%f\t%d\n",
      head->population_id,
      seq_region_id,
      head->position,
      head->var_id,
      next->position,
      next->var_id,
      stats.r2,
      fabs(stats.d_prime),
      stats.N
    );
  }
}

void usage(char *prog) {
  fprintf(stderr, "Usage: %s -f [input.vcf.gz] -r [chr:start-end] -l [optional_sample_list] (-g [input_two.vcf.gz] -s [chr:start-end]) > output.txt\n", prog);
}

int main(int argc, char *argv[]) {

  // parse args
  int c;
  char *files[2];
  char *regions[2];
  char *samples_list;
  char *variant = "";
  int numfiles = 0;
  int numregions = 0;
  int windowsize = WINDOW_SIZE;

  while(1) {
    static struct option long_options[] = {
      {"file",    required_argument, 0, 'f'},
      {"region",  required_argument, 0, 'r'},
      {"file2",   required_argument, 0, 'g'},
      {"region2", required_argument, 0, 's'},
      {"samples", required_argument, 0, 'l'},
      {"window",  required_argument, 0, 'w'},
      {"variant", required_argument, 0, 'v'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "f:g:l:r:s:w:v:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {

      case 'f':
        files[numfiles++] = optarg;
        break;

      case 'g':
        files[numfiles++] = optarg;
        break;

      case 'r':
        regions[numregions++] = optarg;
        break;

      case 's':
        regions[numregions++] = optarg;
        break;

      case 'l':
        samples_list = optarg;
        break;

      case 'w':
        windowsize = atoi(optarg);

      case 'v':
        variant = optarg;
        break;

      case '?':
        /* getopt_long already printed an error message. */
        break;

      default:
        abort ();
    }
  }

  if(numfiles == 0) {
    fprintf(stderr, "No file(s) specified with -f/-g\n");
    usage(argv[0]);
    return USER_ERROR;
  }

  if(numregions == 0) {
    fprintf(stderr, "No region(s) specified with -r/-s\n");
    usage(argv[0]);
    return USER_ERROR;
  }

  if(numfiles != numregions) {
    fprintf(stderr, "Number of files does not match number of regions\n");
    usage(argv[0]);
    return USER_ERROR;
  }
  if(numfiles > 1) {
    windowsize = 1000000000;
  }

  // init vars
  FILE *fh;

  int position;
  int personid;
  int seq_region_id;
  int population_id = 1;
  char genotype[2];

  Locus_list locus_list;
  Locus_info* l_tmp;

  // open output
  fh = stdout; // fopen("output.txt","w");

  init_locus_list(&locus_list);

  int f;
  for(f=0; f<numfiles; f++) {

    // open htsFile and index
    htsFile *htsfile = hts_open(files[f], "rz");
    tbx_t *idx = tbx_index_load(files[f]);

    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(htsfile);

    // use sample list if provided
    // this speeds up VCF parsing
    if(samples_list) {

      // can be a file
      int is_file = 1;

      // or a comma-separated list
      if(strstr(samples_list, ",") != NULL) {
        is_file = 0;
      }

      if(bcf_hdr_set_samples(hdr, samples_list, is_file) < 0) {
        fprintf(stderr, "Failed to read or set samples\n");
        return USER_ERROR;
      }
    }

    // query
    hts_itr_t *iter = tbx_itr_querys(idx, regions[f]);

    // dive out without iter
    if(!iter) return 0;

    // set up vars
    kstring_t str = {0,0,0};
    bcf1_t *line = bcf_init();
    int ngt, *gt_arr = NULL, ngt_arr = 0, i = 0, j = 0, allele;

    // iterate over file
    vcf_line: while(tbx_itr_next(htsfile, idx, iter, &str) > 0) {

      // parse into vcf struct as line
      if(vcf_parse(&str, hdr, line) == 0) {

        seq_region_id = (int) line->rid;

        // get position
        // in htslib it's 0-indexed
        // we also want to increment by one if it's not a SNP
        // to account for the base added to REF and ALT in VCF format
        position = line->pos + (2 - bcf_is_snp(line));

        // get variant id
        // have to do a string copy otherwise the reference to the last one gets passed around indefinitely
        bcf_unpack(line, 1);
        bcf_dec_t *d = &line->d;      
        char *var_id = malloc(strlen(d->id)+1);
        strcpy(var_id, d->id);

        // get genotypes
        ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);

        if(ngt > 0) {

          // quickly scan through for non-ref alleles
          // this way we can exclude non-variant sites before doing any analysis
          int has_alt = 0;
          for(i=0; i<ngt; i++) {
            if(gt_arr[i] > 3) {
              has_alt = 1;
              break;
            }
          }
          if(has_alt == 0) goto vcf_line;

          // gt_arr is an array of alleles
          // we need to break this down into per-sample sets to turn into genotypes
          int alleles_per_gt = ngt / line->n_sample;

          // for now skip unless ploidy == 2
          if(alleles_per_gt != 2) goto vcf_line;

          // iterate over genotypes
          for(i=0; i<line->n_sample; i++) {
            genotype: personid = i + 1;

            // iterate over alleles
            for(j=0; j<alleles_per_gt; j++) {

              // convert hts encoding to a plain int
              // 0 = missing
              // 1 = REF
              // 2 = ALT1
              // 3 = ALT2 etc
              allele = gt_arr[(alleles_per_gt*i) + j]/2;

              // for now we'll only deal with REF or ALT1
              if(allele == 1) {
                genotype[j] = 'A';
              }
              else if(allele == 2) {
                genotype[j] = 'a';
              }
              // if any alleles are missing or ALT2+ just skip this whole variant
              else {
                i++;
                if(i >= line->n_sample) goto vcf_line;
                goto genotype;
              }
            }

            /* Check both are 'a' or 'A' */
            if ((genotype[0] | genotype[1] | 0x20) != 'a') {
              fprintf(stderr, "Genotype must be AA, Aa or aa, not %d\t%d\n", position, personid);
              return USER_ERROR;
            }

            /* Make all hets the same order */
            if (genotype[0] == 'a' && genotype[1] == 'A') {
              genotype[0] = 'A'; genotype[1] = 'a';
            }

            // init locus using enqueue on first genotype
            if(i == 0) {
              enqueue(&locus_list, position, var_id, population_id, personid, genotype2int(genotype));

              // get l_tmp ref to use for subsequent genotypes
              l_tmp = &locus_list.locus[locus_list.tail];
            }

            // subsequent genotypes get added to l_tmp
            else {
              l_tmp->genotypes[l_tmp->number_genotypes].person_id = personid;
              l_tmp->genotypes[l_tmp->number_genotypes].genotype = genotype2int(genotype);
              l_tmp->number_genotypes++;
              if (l_tmp->number_genotypes == SIZE) {
                fprintf(stderr, "Number of genotypes supported by the program (%d) exceeded\n", SIZE);
                return SYSTEM_ERROR;
              }
            }
          }
        }

        /*check if the new position is farther than the limit.*/
        /*if so, calculate the ld information for the values in the array*/
        while(
          (locus_list.tail >= locus_list.head) &&
          (abs(locus_list.locus[locus_list.head].position - position) > windowsize)
        ) {
//          fprintf(stderr, "Calculate LD 1\n");

          calculate_ld(&locus_list, seq_region_id, fh, variant);
          dequeue(locus_list);  
        }
        if (locus_list.tail < locus_list.head) {
          /* Can reset the queue to the beginning */
          locus_list.head = 0;
          locus_list.tail = -1;
        }
      }
    }

    tbx_destroy(idx);
    tbx_itr_destroy(iter);
    bcf_hdr_destroy(hdr);
    hts_close(htsfile);
  }

  l_tmp = &locus_list.locus[locus_list.tail];

  while(locus_list.tail >= locus_list.head){

    calculate_ld(&locus_list, seq_region_id, fh, variant);
    dequeue(locus_list);
  }

  return 0;
}

