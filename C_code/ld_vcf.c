/* A program for calculating linkage disequilibrium stats */

/*
 Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
#include <unistd.h>

#include "tbx.h"
#include "vcf.h"

#define SIZE 1000000
#define WINDOW_SIZE 100000
#define INITIAL_LIST_SIZE 256
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

void dequeue(Locus_list *ll) {
  ll->head++;
}

void enqueue(Locus_list *ll, int pos, char *var_id, int personid, uint8_t genotype) {
  Locus_info *l;

  ll->tail++;
  if (ll->tail == ll->sz)
    reallocate_locus_list(ll);

  l = &ll->locus[ll->tail];
  l->position = pos;
  l->var_id = var_id;
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

  // fprintf(
  //   stderr,
  //   "AABB %i\nAABb %i\nAAbb %i\nAaBB %i\nAaBb %i\nAabb %i\naaBB %i\naaBb %i\naabb %i\n",
  //   AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb
  // );

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
    // fprintf(stderr, "theta %f thetaprev %f\n", theta, thetaprev);
  }

  /*now calculate stats*/
  major_freqs(&haplotypes,&f_A,&f_B);
  // fprintf(stderr, "f_A %f\tf_B %f\n", f_A, f_B);
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

void calculate_ld(const Locus_list *ll, FILE *fh, char *variant, int windowsize){
  Locus_info *next, *head;
  Stats stats;
  
  /* Doesn't look like it, but sets head and next to the first and
     second entries - I love C, sometimes */

  next = &ll->locus[ll->head];
  head = next++;
  int i;
  for (i = ll->head; i < ll->tail; i++, next++) {

    // skip if > windowsize
    if(windowsize > 0 && abs(head->position - next->position) > windowsize)
      continue;

    // skip if user has provided variant and neither matches ID
    if(variant[0] != 0 && strcmp(head->var_id, variant) != 0 && strcmp(next->var_id, variant) != 0)
      continue;

    calculate_pairwise_stats(head, next, &stats);
    if ((float) stats.r2 < 0.05 || stats.N < 40 || (float) stats.r2 > 1 || (float) stats.d_prime > 1)
      continue;
    fprintf(fh, "%d\t%d\t%d\t%s\t%d\t%s\t%f\t%f\t%d\n",
      1,   // this used to be population_id, but we're ignoring that now. Placeholder for file format compatiblity
      1,   // this used to be seq_region_id, but we're ignoring that now. Placeholder for file format compatiblity
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

int get_genotypes(Locus_list *locus_list, char *variant, bcf_hdr_t *hdr, bcf1_t *line) {

  Locus_info* l_tmp;
  int ngt, *gt_arr = NULL, ngt_arr = 0, i = 0, j = 0, allele, personid;

  // get position
  // in htslib it's 0-indexed
  // we also want to increment by one if it's not a SNP
  // to account for the base added to REF and ALT in VCF format
  int position = line->pos + (2 - bcf_is_snp(line));

  // get variant id
  // have to do a string copy otherwise the reference to the last one gets passed around indefinitely
  bcf_unpack(line, 1);
  bcf_dec_t *d = &line->d;      
  char *var_id = malloc(strlen(d->id)+1);
  strcpy(var_id, d->id);

  // get genotypes
  ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);

  if(ngt == 0) return 0;

  // quickly scan through for non-ref alleles
  // this way we can exclude non-variant sites before doing any analysis
  int has_alt = 0;
  for(i=0; i<ngt; i++) {
    if(gt_arr[i] > 3) {
      has_alt = 1;
      break;
    }
  }
  if(has_alt == 0) return 0;

  // gt_arr is an array of alleles
  // we need to break this down into per-sample sets to turn into genotypes
  int alleles_per_gt = ngt / line->n_sample;

  // for now skip unless ploidy == 2
  if(alleles_per_gt != 2) return 0;

  int initialised = 0;

  // iterate over genotypes
  for(i=0; i<line->n_sample; i++) {
    genotype: personid = i + 1;

    char genotype[2];

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
        if(i >= line->n_sample) return 0;
        goto genotype;
      }
    }

    /* Check both are 'a' or 'A' */
    if ((genotype[0] | genotype[1] | 0x20) != 'a') {
      fprintf(stderr, "Genotype must be AA, Aa or aa, not %d\t%d\n", position, personid);
      return EXIT_FAILURE;
    }

    /* Make all hets the same order */
    if (genotype[0] == 'a' && genotype[1] == 'A') {
      genotype[0] = 'A'; genotype[1] = 'a';
    }

    // init locus using enqueue on first genotype
    if(initialised == 0) {
      enqueue(locus_list, position, var_id, personid, genotype2int(genotype));

      // get l_tmp ref to use for subsequent genotypes
      l_tmp = &locus_list->locus[locus_list->tail];

      initialised = 1;
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

  return position;
}

void process_window(Locus_list *locus_list, int windowsize, char *variant, FILE *fh, int position) {
  /*check if the new position is farther than the limit.*/
  /*if so, calculate the ld information for the values in the array*/
  while(
    (locus_list->tail >= locus_list->head) &&
    (abs(locus_list->locus[locus_list->head].position - position) > windowsize)
  ) {
    calculate_ld(locus_list, fh, variant, windowsize);
    dequeue(locus_list);  
  }
  if (locus_list->tail < locus_list->head) {
    /* Can reset the queue to the beginning */
    locus_list->head = 0;
    locus_list->tail = -1;
  }
}

void usage(char *prog) {
  fprintf(stderr, "Usage: %s -f [input.vcf.gz] -r [chr:start-end] -l [optional_sample_list] (-g [input_two.vcf.gz] -s [chr:start-end]) > output.txt\n", prog);
}

char** read_variants_file(char *variants_file) {
  int lines_allocated = 128;
  int max_line_len = 100;
  char** include_variants = (char **)malloc(sizeof(char*)*lines_allocated);
  if (include_variants==NULL) {
    fprintf(stderr, "Out of memory.\n");
    exit(SYSTEM_ERROR);
  }

  FILE *in;
  if ((in = fopen(variants_file, "r"))==NULL) {
    perror("Could not open input file");
    exit(SYSTEM_ERROR);
  }

  int l;
  for(l=0; 1; l++) {
    int j;

    // increase memory as required
    if (l >= lines_allocated) {
      int new_size;

      /* Double our allocation and re-allocate */
      new_size = lines_allocated*2;
      include_variants = (char **)realloc(include_variants,sizeof(char*)*new_size);
      if (include_variants==NULL) {
        fprintf(stderr,"Out of memory.\n");
        exit(SYSTEM_ERROR);
      }
      lines_allocated = new_size;
    }

    include_variants[l] = malloc(max_line_len);

    if(include_variants[l]==NULL) {
      fprintf(stderr,"Out of memory.\n");
      exit(SYSTEM_ERROR);
    }
    if(fgets(include_variants[l], max_line_len-1, in)==NULL) break;

    /* Get rid of CR or LF at end of line */
    for(
      j=strlen(include_variants[l])-1;
      j>=0 && (include_variants[l][j]=='\n' || include_variants[l][j]=='\r');
      j--
    ) { 
      include_variants[l][j]='\0';
    }
  }

  fclose(in);

  return include_variants;
}

int check_include_variants(bcf1_t *line, char** include_variants, char* variant) {
  bcf_unpack(line, 1);
  bcf_dec_t *d = &line->d;      
  char *var_id = malloc(strlen(d->id)+1);
  strcpy(var_id, d->id);

  // could be the variant given with -v
  if((variant[0] != '\0') && strcmp(var_id, variant) == 0) return 1;

  // or could be in the file given
  int i;
  for(i=0; include_variants[i][0] != '\0'; i++) {
    if(strcmp(include_variants[i], var_id) == 0) return 1;
  }

  return 0;
}

int main(int argc, char *argv[]) {

  // parse args
  int c;
  char *files[2];
  char *regions[2];
  char *samples_list;
  char *variants_file;
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
      {"include_variants", required_argument, 0, 'n'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "f:g:l:r:s:w:v:n:", long_options, &option_index);

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
        windowsize = (int) atoi(optarg);
        break;

      case 'v':
        variant = optarg;
        break;

      case 'n':
        variants_file = optarg;
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
    return EXIT_FAILURE;
  }

  if(numregions == 0) {
    fprintf(stderr, "No region(s) specified with -r/-s\n");
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  if(numfiles != numregions) {
    fprintf(stderr, "Number of files does not match number of regions\n");
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  if(numfiles > 1) {
    windowsize = 1000000000;
  }

  // variant list in file
  char **include_variants;
  int have_include_variants = 0;
  // include_variants[0][0] = '\0';
  if(access( variants_file, F_OK) != -1 ) {
    include_variants = read_variants_file(variants_file);
    have_include_variants = 1;
  }

  // open output
  FILE *fh;
  fh = stdout; // fopen("output.txt","w");

  // init vars
  Locus_list locus_list;
  init_locus_list(&locus_list);
  int f, position;

  for(f=0; f<numfiles; f++) {

    // open htsFile
    htsFile *htsfile = hts_open(files[f], "rz");

    if(!htsfile) {
      fprintf(stderr, "Unable to open file %s\n", files[f]);
      return EXIT_FAILURE;
    }

    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(htsfile);

    if(!hdr) {
      fprintf(stderr, "Unable to read header from file %s\n", files[f]);
      return EXIT_FAILURE;
    }

    // use sample list if provided
    // this speeds up VCF parsing
    if(samples_list) {

      // can be a file
      int is_file = 1;

      // or a comma-separated list
      if(strstr(samples_list, ",") != NULL) {
        is_file = 0;
      }

      if(access( samples_list, F_OK ) != -1  && bcf_hdr_set_samples(hdr, samples_list, is_file) < 0) {
        fprintf(stderr, "Failed to read or set samples\n");
        return EXIT_FAILURE;
      }
    }

    // get file format and act accordingly
    enum htsExactFormat format = hts_get_format(htsfile)->format;

    if(format == vcf) {
    
      // open index
      tbx_t *idx = tbx_index_load(files[f]);

      if(!idx) {
        fprintf(stderr, "Could not load .tbi/.csi index for file %s\n", files[f]);
        return EXIT_FAILURE;
      }

      // query
      hts_itr_t *itr = tbx_itr_querys(idx, regions[f]);

      // dive out without iter
      if(!itr) return 0;

      // set up vars
      kstring_t str = {0,0,0};
      bcf1_t *line = bcf_init();

      // iterate over file
      while(tbx_itr_next(htsfile, idx, itr, &str) > 0) {

        // parse into vcf struct as line
        if(vcf_parse(&str, hdr, line) == 0) {

          // check include_variants
          if(have_include_variants && check_include_variants(line, include_variants, variant) == 0) continue;

          position = get_genotypes(&locus_list, variant, hdr, line);
          if(position > 0) process_window(&locus_list, windowsize, variant, fh, position);
        }
      }

      free(str.s);
      tbx_destroy(idx);
      tbx_itr_destroy(itr);
    }

    else if(format == bcf) {

      // open index
      hts_idx_t *idx = bcf_index_load(files[f]);

      if(!idx) {
        fprintf(stderr, "Could not load .csi index for file %s\n", files[f]);
        return EXIT_FAILURE;
      }

      // query
      hts_itr_t *itr = bcf_itr_querys(idx, hdr, regions[f]);

      // dive out without iter
      if(!itr) return 0;

      bcf1_t *line = bcf_init();

      while(bcf_itr_next(htsfile, itr, line) >= 0) {

        // check include_variants
        if(have_include_variants && check_include_variants(line, include_variants, variant) == 0) continue;

        position = get_genotypes(&locus_list, variant, hdr, line);
        if(position > 0) process_window(&locus_list, windowsize, variant, fh, position);
      }

      hts_idx_destroy(idx);
      bcf_itr_destroy(itr);
    }

    else {
      fprintf(stderr, "Unsupported format for file %s\n", files[f]);
      return EXIT_FAILURE;
    }

    bcf_hdr_destroy(hdr);

    if ( hts_close(htsfile) ) {
      fprintf(stderr, "hts_close returned non-zero status: %s\n", files[f]);
      return EXIT_FAILURE;
    }
  }

  // process any remaining buffer
  process_window(&locus_list, 0, variant, fh, position);

  return 0;
}

