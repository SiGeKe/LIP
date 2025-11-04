// Header-Files

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

// Definitions

#define N 10001 // max. number of lattice sites
#define M 21 // max. number of lipid types
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}

// Variables

int reverse[4];

double temperature;
int lattice_size, length[4];
int nr_lattice;
int cholpos[N];
int particlet[N], order[N];
int cmax_dd,cmax_cd,cmax_cc;

int neigh[N][7], neigh_cd[N][7], neigh_dc[N][7];
int neigh2_dd[N][4][10],neigh2_cd[N][4][10],neigh2_cc[N][4][10];

int t_total;
int nr_data, nr_count, nr_average;
double pequil; 

int nr_types;
double c_chol, c_lip[M];
int nr_chol,nrt[M];

double H_self[M][151], H_dd[M][M][151][10], H_cd[M][151][10], H_cc[10];
double Entro[M][151];
double nneigh_dd[M][5], nneigh_cd[M][5];

double H_dd_avg[M][M][301][10]; // ACHTUNG

int neightype[M][M], neightype_dc[M], neightype_cc;
double p_MC[M][151],p_MC_neigh[M][151][5];

double p_MD[151];

double init_origin[M][151];

char project_name[1024];
double read_temperature;
int mv_type, inter_type;
int dlipc;
int orderwidth,dim;

double quotf;

int first_runs;

double a_S;
double b_S;
double c_N;

int cc_contacts;

double N_c_mean;

// Functions

void fisher_yate_shuffle(int *);
int lattice_sites[N+1];

void main();

void read_input();

void calc_all();
void calc_evolution();
void make_sign();

void make_neighbors();
void make_particle_numbers();
void make_initialization();
void make_interaction();
void lattice_init();
void sort(int,double *,int *);

void mc_move_order(int);
void mc_move_particle(int);
void mc_move_chol(int);

void plot_info();
void plot_order(int,char[]);
void plot_order_neigh(int,char[]);
void plot_lattice(int,char[]);
void plot_Entro(int,char[]);

void count_neighbors();
int count_cc(int);
int count_dc(int);
int count2_dd(int,int);
int count2_dc(int,int);
int count2_cc(int,int);
void count_order();
void norm_order();

void norm_Entro();

int min_order(int,int);
int max_order(int,int);
int avg_order(int,int);
double avg_part(int);
int per(int,int);

void MD_order();
void update_Entro(double);
double order_diff();
double order_diff_neigh(int);
double P_init();
double L_init();

double ener_dd(int,int);
double ener_cd(int,int);
double ener_cc(int,int);
double ener_self(int);

double NN_CD_L(int);
double NN_CD_P(int);
double NN_DD_L(int);
double NN_DD_P(int);

double H_self_L(int);
double H_self_P(int);

double H_DD_LL(int,int);
double H_DD_PP(int,int);
double H_DD_LP(int,int);

double H_DD_PP_avg(int,int); 
double H_DD_LL_avg(int,int);

double H_CC_L(int);
double H_CC_P(int);

double H_CD_L(int,int);
double H_CD_P(int,int);

double E_L(int);
double E_P(int);

double norm_dist(double,double,double);

double avg_order_mc();
double avg_order_md();
int max_order_mc();
int max_order_md();
double avg_order_neigh(int);

// Main Functions

void main()
{
    int seed = 128;

    read_input(); 

    reverse[0] = 1, reverse[1] = 0, reverse[2] = 3, reverse[3] = 2;

    MD_order();

    if (nr_types == 1) c_lip[1] = 1;
    c_chol = (c_chol)/(1-c_chol);
    temperature = read_temperature/120.0;
    length[1] = length[2] = lattice_size;

    inter_type = 1; 
    nr_data = 10; 

    srandom(seed);
    quotf = pow(2.0,-31.0);

    calc_all();

    if (dlipc == 1 && nr_types == 1) printf("Binary System - DLiPC\n");
    if (dlipc != 1 && nr_types == 1) printf("Binary System - DPPC\n");
    if (nr_types == 2) printf("Ternary System - DLiPC/DPPC\n");
    printf("\n");
    printf("Chol.-Conc. on lattice: %lf\n", c_chol);
    c_chol = (c_chol)/(1+c_chol);
    printf("Chol.-Conc. in system: %lf\n", c_chol);
    printf("Move-Type: %d\n", mv_type);
    printf("Temperature: %lf\n", read_temperature);

    plot_info();
}

void calc_all()
{
    int run = 1;
    char ufname[1024];

    double diff_10 = 0.0;
    double diff_20 = 0.0;
    double diff_30 = 0.0;

    double avg_mc_10 = 0.0;
    double avg_mc_20 = 0.0;
    double avg_mc_30 = 0.0;

    double avg_md_10 = 0.0;
    double avg_md_20 = 0.0;
    double avg_md_30 = 0.0;

    double max_md_10 = 0.0;
    double max_mc_10 = 0.0;

    double avg_10_0 = 0.0;
    double avg_10_1 = 0.0;
    double avg_10_4 = 0.0;

    double avg_20_0 = 0.0;
    double avg_20_1 = 0.0;
    double avg_20_4 = 0.0;

    double avg_30_0 = 0.0;
    double avg_30_1 = 0.0;
    double avg_30_4 = 0.0;

    double cc_10 = 0.0;
    double cc_20 = 0.0;
    double cc_30 = 0.0;

    double diff_10_1 = 0.0;
    double diff_10_4 = 0.0;
    
    double diff_20_1 = 0.0;
    double diff_20_4 = 0.0;

    double diff_30_1 = 0.0;
    double diff_30_4 = 0.0;

    double N_c_mean_10 = 0.0;
    double N_c_mean_20 = 0.0;
    double N_c_mean_30 = 0.0;

    a_S = 0.0;
    b_S = 0.0;
    c_N = 0.0;

    double alpha = 0.0;

    while (0==0){

        /*---10%---*/
        c_chol = 0.1;
        MD_order();
        c_chol = (c_chol)/(1-c_chol);

        make_neighbors();
        make_particle_numbers();
        make_initialization();
        make_interaction();

        calc_evolution();

        diff_10 = order_diff()*10000;
        diff_10_1 = order_diff_neigh(1)*10000;
        diff_10_4 = order_diff_neigh(4)*10000;
        avg_mc_10 = avg_order_mc();
        avg_md_10 = avg_order_md();
        max_mc_10 = max_order_mc();
        max_md_10 = max_order_md();
        avg_10_0 = avg_order_neigh(0);
        avg_10_1 = avg_order_neigh(1);
        avg_10_4 = avg_order_neigh(4);
        count_neighbors();
        cc_10 = (double)(cc_contacts)/((double)t_total - pequil*(double)t_total);

        N_c_mean_10 = N_c_mean;

        plot_order(run,ufname);

        /*---20%---*/
        c_chol = 0.2;
        MD_order();
        c_chol = (c_chol)/(1-c_chol);

        make_neighbors();
        make_particle_numbers();
        make_initialization();
        make_interaction();

        calc_evolution();

        diff_20 = order_diff()*10000;
        diff_20_1 = order_diff_neigh(1)*10000;
        diff_20_4 = order_diff_neigh(4)*10000;
        avg_mc_20 = avg_order_mc();
        avg_md_20 = avg_order_md();
        avg_20_0 = avg_order_neigh(0);
        avg_20_1 = avg_order_neigh(1);
        avg_20_4 = avg_order_neigh(4);
        count_neighbors();
        cc_20 = (double)(cc_contacts)/((double)t_total - pequil*(double)t_total);

        N_c_mean_20 = N_c_mean;

        plot_order(run,ufname);

        /*---30%---*/
        c_chol = 0.3;
        MD_order();
        c_chol = (c_chol)/(1-c_chol);

        make_neighbors();
        make_particle_numbers();
        make_initialization();
        make_interaction();

        calc_evolution();

        diff_30 = order_diff()*10000;
        diff_30_1 = order_diff_neigh(1)*10000;
        diff_30_4 = order_diff_neigh(4)*10000;
        avg_mc_30 = avg_order_mc();
        avg_md_30 = avg_order_md();
        avg_30_0 = avg_order_neigh(0);
        avg_30_1 = avg_order_neigh(1);
        avg_30_4 = avg_order_neigh(4);
        count_neighbors();
        cc_30 = (double)(cc_contacts)/((double)t_total - pequil*(double)t_total);

        N_c_mean_30 = N_c_mean;

        plot_order(run,ufname);

        printf("Run: %d\n", run);
        printf("10 - Difference: %lf, MD: %lf, MC: %lf, Contacts: %lf, <N_c>: %lf\n", diff_10, avg_md_10, avg_mc_10, cc_10, N_c_mean_10);
        printf("20 - Difference: %lf, MD: %lf, MC: %lf, Contacts: %lf, <N_c>: %lf\n", diff_20, avg_md_20, avg_mc_20, cc_20, N_c_mean_20);
        printf("30 - Difference: %lf, MD: %lf, MC: %lf, Contacts: %lf, <N_c>: %lf\n", diff_30, avg_md_30, avg_mc_30, cc_30, N_c_mean_30);
        printf("All - Difference: %lf, a_S: %lf, b_S: %lf, c_N: %lf, alpha: %lf", (diff_10+diff_20+diff_30/3.0),a_S,b_S,c_N,alpha);
        printf("\n\n");

        /*Iterations*/

        a_S -= 0.1;
        b_S -= 0.5;
        c_N -= 0.01;

        alpha *= 0.9;


        

        run ++;

        if (run == 10) break;
    }
}

void calc_evolution()
{
    FILE *fconfig;
    int l,n,t;
    char fname[1024], str[12];
    char ufname[1024];

    for (n=0;n<=nr_average;n++){

        lattice_init();
        make_initialization();
        count_neighbors();

        strcpy(fname, project_name);
        strcat(fname, "_contacts_");
        sprintf(str, "%d", n);
        strcat(fname,str);
        strcat(fname,".dat");

        for (t=1;t<=t_total;t++){

            if ((double)t > pequil*(double)t_total){
                
                count_order();
                count_neighbors();
                cc_contacts += neightype_cc;

            }

            fisher_yate_shuffle(lattice_sites);

            double N_c = 0.0;

            for (l=0;l<nr_lattice;l++){
                int lm = lattice_sites[l];
                mc_move_order(lm);
                mc_move_particle(lm);
                if (cholpos[lm] == 1) mc_move_chol(lm);

                N_c += (double)count_dc(lm);
            }
            N_c_mean += N_c/(double)nr_lattice;
        }
        N_c_mean /= (double)t_total;

        norm_order();

    }
}

void make_particle_numbers()
{
    int i;
    int sum = 0;
    double num = 0.0;

    for (i=1;i<nr_types;i++){
        num = (double)(nr_lattice)*c_lip[i];
        nrt[i] = (int)round(num);
        sum += nrt[i];
    }
    nrt[nr_types] = nr_lattice - sum;
}

void make_initialization()
{
    int i,s,c;
    for (i=1;i<=nr_types;i++){
        for (s=0;s<=150;s++){
            p_MC[i][s] = 0;
            for (c=0;c<=4;c++){
                p_MC_neigh[i][s][c] = 0;
            }
        }
    }

    cc_contacts = 0.0;
    N_c_mean = 0.0;

    for (int l=0; l<=nr_lattice;l++){
        lattice_sites[l] = l;
    }
}

void make_interaction()
{
    int a,b,c;
    int s;

    if (inter_type == 1){
        if (dlipc == 1 && nr_types == 1){
            for (a=0;a<=4;a++){
                nneigh_dd[1][a] = NN_DD_L(a);
                nneigh_cd[1][a] = NN_CD_L(a);
            }
            for (s=0;s<=150;s++){
                H_self[1][s] = H_self_L(s);
                Entro[1][s] = E_L(s);
                for (b=0;b<=6;b++){
                    H_dd[1][1][s][b] = H_DD_LL(s,b);
                    H_cc[b] = H_CC_L(b);
                }
                for (c=0;c<=5;c++){
                    H_cd[1][s][c] = H_CD_L(s,c);
                }
            }
            for (s=0;s<=300;s++){
                for (b=0;b<=6;b++){
                    H_dd_avg[1][1][s][b] = H_DD_LL_avg(s,b);
                }
            }
        }
        if (dlipc != 1 && nr_types == 1){
            for (a=0;a<=4;a++){
                nneigh_dd[1][a] = NN_DD_P(a);
                nneigh_cd[1][a] = NN_CD_P(a);
            }
            for (s=0;s<=150;s++){
                H_self[1][s] = H_self_P(s);
                Entro[1][s] = E_P(s);
                for (b=0;b<=6;b++){
                    H_dd[1][1][s][b] = H_DD_PP(s,b); 
                    H_cc[b] = H_CC_P(b);
                }
                for (c=0;c<=5;c++){
                    H_cd[1][s][c] = H_CD_P(s,c);
                }
            }
            for (s=0;s<=300;s++){
                for (b=0;b<=6;b++) {
                    H_dd_avg[1][1][s][b] = H_DD_PP_avg(s,b);
                }
            }
        }
        if (dlipc == 1 && nr_types == 2){
            for (a=0;a<=4;a++){
                nneigh_dd[1][a] = NN_DD_L(a);
                nneigh_dd[2][a] = NN_DD_P(a);
                nneigh_cd[1][a] = NN_CD_L(a);
                nneigh_cd[2][a] = NN_CD_P(a);
            }
            for (s=0;s<=150;s++){
                H_self[1][s] = H_self_L(s);
                H_self[2][s] = H_self_P(s);
                Entro[1][s] = E_L(s);
                Entro[2][s] = E_P(s);
                for (b=0;b<=6;b++){
                    H_dd[1][1][s][b] = H_DD_LL(s,b);
                    H_dd[2][2][s][b] = H_DD_PP(s,b);
                    H_dd[1][2][s][b] = H_dd[2][1][s][b] = H_DD_LP(s,b);
                    H_cc[b] = H_CC_P(b);
                }
                for (c=0;c<=5;c++){
                    H_cd[1][s][c] = H_CD_L(s,c);
                    H_cd[1][s][c] = H_CD_P(s,c);
                }
            }
        }
    }
}


void lattice_init()
{
  int sum,i,j,k;
  double x[N];
  for (i=1 ; i <= nr_lattice ; i++)
  x[i] = random();

  i = 1;
  for (j=1 ; j <= nr_types ; j++)
  for (k=1 ; k <= nrt[j] ; k++)
  {
    particlet[i] = j;
    i++;
  }
  sort(nr_lattice,x,particlet);
  for (i=0 ; i < nr_lattice ; i++)
  particlet[i] = particlet[i+1];
  for (i=0 ; i < nr_lattice; i++) 
  {
    if (dlipc == 1 && particlet[i] == 1){
        while (order[i] > 120 || order[i] < 30){ 
          order[i] = 50 + (int)(30*random()*quotf);
        }  
    }
    else{
      while (order[i] > 150 || order[i] < 30){
      order[i] = 30 + (int)(30*(random()*quotf));
      }
    } 
  }

  for (i=1 ; i <= nr_chol; i++)
  {
    cholpos[i] = 1; x[i] = random(); 
  }
  for (i=nr_chol+1 ; i <= nr_lattice; i++)
  {
    cholpos[i] = 0; x[i] = random(); 
  }
  sort(nr_lattice,x,cholpos);
  for (i=0 ; i < nr_lattice ; i++)
  cholpos[i] = cholpos[i+1];
}

// MC-Moves

void mc_move_order(int l)
{
    double delta,jc_old,jc_new;
    int delta_order;
    int m;

    jc_old = 0;
    for (m=0;m<2*dim;m++){
        jc_old += ener_dd(l,m);
        if (cholpos[neigh_dc[l][m]] == 1){
            jc_old += ener_cd(neigh_dc[l][m],reverse[m]);
        }
    }
    jc_old += ener_self(l);

    delta_order = random()%orderwidth - (orderwidth-1)/2;
    order[l] += delta_order;
    if (dlipc == 1 && particlet[l] == 1){
        if (order[l] < 30 || order[l] > 120){
            order[l] -= delta_order;
            return;
        }
    }
    else{
        if (order[l] < 30 || order[l] > 150){
            order[l] -= delta_order;
            return;
        }
    }


    jc_new = 0;
    for (m=0;m<2*dim;m++){
        jc_new += ener_dd(l,m);
        if (cholpos[neigh_dc[l][m]] == 1){
            jc_new += ener_cd(neigh_dc[l][m],reverse[m]);
        }
    }
    jc_new += ener_self(l);

    delta = jc_new - jc_old;
    if (exp(-delta/temperature) < random()*quotf){
        order[l] -= delta_order;
    }
}

void mc_move_particle(int l)
{
    double delta,jc_old,jc_new;
    int m,l_new;

    if (mv_type == 0){
        l_new = random()%(nr_lattice);
        if (l_new == l) l_new = (l_new + 1)%nr_lattice;
    }
    if (mv_type == 1){
        m = random()%(2*dim);
        l_new = neigh[l][m];
    }

    if (particlet[l] != particlet[l_new]){

        jc_old = 0;
        for (m=0;m<2*dim;m++){
            if (neigh[l][m] != l_new) jc_old += ener_dd(l,m);
            else jc_old += 0.5*ener_dd(l,m);

            if (neigh[l_new][m] != l) jc_old += ener_dd(l_new,m);
            else jc_old += 0.5*ener_dd(l_new,m);

            if (cholpos[neigh_dc[l_new][m]] == 1){
                jc_old += ener_cd(neigh_dc[l_new][m],reverse[m]);
            }
            if (cholpos[neigh_dc[l][m]] == 1){
                jc_old += ener_cd(neigh_dc[l][m],reverse[m]);
            }
        }

        SWAP(particlet[l],particlet[l_new]);
        SWAP(order[l],order[l_new]);

        jc_new = 0;
        for (m=0;m<2*dim;m++){
            if (neigh[l][m] != l_new) jc_new += ener_dd(l,m);
            else jc_new += 0.5*ener_dd(l,m);

            if (neigh[l_new][m] != l) jc_new += ener_dd(l_new,m);
            else jc_new += 0.5*ener_dd(l_new,m);

            if (cholpos[neigh_dc[l_new][m]] == 1){
                jc_new += ener_cd(neigh_dc[l_new][m],reverse[m]);
            }
            if (cholpos[neigh_dc[l][m]] == 1){
                jc_new += ener_cd(neigh_dc[l][m],reverse[m]);
            }
        }

        delta = jc_new - jc_old;

        if (exp(-delta/temperature) < random()*quotf){
            SWAP(particlet[l],particlet[l_new]);
            SWAP(order[l],order[l_new]);
        }
    }
}

void mc_move_chol(int l)
{
    int o,k;
    int m;
    int x,y,g,h;
    double delta,jc_old,jc_new;
    int l_new;

    int inter_cc[nr_lattice][4],inter_cd[nr_lattice][4],inter_dd[nr_lattice][4];
    int inter_self[nr_lattice];

    if (mv_type == 0){
        l_new = random()%(nr_lattice);
        if (l_new == l) l_new = (l_new + 1)%nr_lattice;
    }
    if (mv_type == 1){
        m = random()%(2*dim);
        l_new = neigh[l][m];
    }

    if (cholpos[l_new] == 0){

        for (o=0;o<nr_lattice;o++){
            inter_self[o] = 0;
            for (m=0;m<4;m++){
                inter_cc[o][m] = 0;
                inter_cd[o][m] = 0;
                inter_dd[o][m] = 0;
            }
        }
        for (m=0;m<2*dim;m++){

            x = neigh[l][m];
            y = neigh[l_new][m];
            g = neigh_cd[l][m];
            h = neigh_cd[l_new][m];

            inter_self[g] ++;   // self_interaction
            inter_self[h] ++;   // self_interaction

            inter_cd[l][m] ++;

            if (cholpos[neigh_dc[g][m]] == 1){
                inter_cd[neigh_dc[g][m]][reverse[m]] ++;
            }
            if (cholpos[neigh_dc[h][m]] == 1){
                inter_cd[neigh_dc[h][m]][reverse[m]] ++;
            }

            for (k=0;k<2*dim;k++){

                if (cholpos[x] == 1){ 
                    inter_cd[x][k] ++;
                    if (cholpos[neigh[x][k]] == 1){
                        inter_cc[x][k] ++;
                        inter_cc[neigh[x][k]][reverse[k]] ++;
                    }
                }

                if (cholpos[y] == 1){
                    inter_cd[y][k] ++;
                    if (cholpos[neigh[y][k]] == 1){
                        inter_cc[y][k] ++;
                        inter_cc[neigh[y][k]][reverse[k]] ++;
                    }
                }

                if (m<2){
                    inter_dd[g][k] ++;
                    inter_dd[neigh[g][k]][reverse[k]] ++;

                    inter_dd[h][k] ++;
                    inter_dd[neigh[h][k]][reverse[k]] ++;
                }
            }

            if (m == 2){
                inter_dd[g][2] ++;
                inter_dd[neigh[g][2]][reverse[2]] ++;
                inter_dd[g][1] ++;
                inter_dd[neigh[g][1]][reverse[1]] ++;

                inter_dd[h][2] ++;
                inter_dd[neigh[h][2]][reverse[2]] ++;
                inter_dd[h][1] ++;
                inter_dd[neigh[h][1]][reverse[1]] ++;
            }

            if (m == 3){
                inter_dd[g][3] ++;
                inter_dd[neigh[g][3]][reverse[3]] ++;
                inter_dd[g][0] ++;
                inter_dd[neigh[g][0]][reverse[0]] ++;

                inter_dd[h][3] ++;
                inter_dd[neigh[h][3]][reverse[3]] ++;
                inter_dd[h][0] ++;
                inter_dd[neigh[h][0]][reverse[0]] ++;
            }
        }

        jc_old = 0;
        for (m=0;m<2*dim;m++){

            x = neigh[l][m];
            y = neigh[l_new][m];
            g = neigh_cd[l][m];
            h = neigh_cd[l_new][m];

            jc_old += ener_self(g);
            jc_old += ener_self(h);

            jc_old += ener_cd(l,m)/(double)(inter_cd[l][m]);

            if (cholpos[neigh_dc[g][m]] == 1){
                jc_old += ener_cd(neigh_dc[g][m],reverse[m])/(double)(inter_cd[neigh_dc[g][m]][reverse[m]]);
            }
            if (cholpos[neigh_dc[h][m]] == 1){
                jc_old += ener_cd(neigh_dc[h][m],reverse[m])/(double)(inter_cd[neigh_dc[h][m]][reverse[m]]);
            }

            for (k=0;k<2*dim;k++){

                if (cholpos[x] == 1){
                    jc_old += ener_cd(x,k)/(double)(inter_cd[x][k]);
                    if (cholpos[neigh[x][k]] == 1) jc_old += ener_cc(x,k)/(double)(inter_cc[x][k]);
                }

                if (cholpos[y] == 1){
                    jc_old += ener_cd(y,k)/(double)(inter_cd[y][k]);
                    if (cholpos[neigh[y][k]] == 1) jc_old += ener_cc(y,k)/(double)(inter_cc[y][k]);
                }

                if (m<2){
                    jc_old += ener_dd(g,k)/(double)(inter_dd[g][k]);
                    jc_old += ener_dd(h,k)/(double)(inter_dd[h][k]);
                }
            }

            if (m == 2){
                jc_old += ener_dd(g,2)/(double)(inter_dd[g][2]);
                jc_old += ener_dd(g,1)/(double)(inter_dd[g][1]);

                jc_old += ener_dd(h,2)/(double)(inter_dd[h][2]);
                jc_old += ener_dd(h,1)/(double)(inter_dd[h][1]);
            }

            if (m == 3){
                jc_old += ener_dd(g,3)/(double)(inter_dd[g][3]);
                jc_old += ener_dd(g,0)/(double)(inter_dd[g][0]);

                jc_old += ener_dd(h,3)/(double)(inter_dd[h][3]);
                jc_old += ener_dd(h,0)/(double)(inter_dd[h][0]);
            }
        }
        
        SWAP(cholpos[l],cholpos[l_new]);

        for (o=0;o<nr_lattice;o++){
            inter_self[o] = 0;
            for (m=0;m<4;m++){
                inter_cc[o][m] = 0;
                inter_cd[o][m] = 0;
                inter_dd[o][m] = 0;
            }
        }
        for (m=0;m<2*dim;m++){

            x = neigh[l][m];
            y = neigh[l_new][m];
            g = neigh_cd[l][m];
            h = neigh_cd[l_new][m];

            inter_self[g] ++;   // ener_self
            inter_self[h] ++;   // ener_self

            inter_cd[l_new][m] ++;

            if (cholpos[neigh_dc[g][m]] == 1){
                inter_cd[neigh_dc[g][m]][reverse[m]] ++;
            }
            if (cholpos[neigh_dc[h][m]] == 1){
                inter_cd[neigh_dc[h][m]][reverse[m]] ++;
            }

            for (k=0;k<2*dim;k++){

                if (cholpos[x] == 1){ 
                    inter_cd[x][k] ++;
                    if (cholpos[neigh[x][k]] == 1){
                        inter_cc[x][k] ++;
                        inter_cc[neigh[x][k]][reverse[k]] ++;
                    }
                }

                if (cholpos[y] == 1){
                    inter_cd[y][k] ++;
                    if (cholpos[neigh[y][k]] == 1){
                        inter_cc[y][k] ++;
                        inter_cc[neigh[y][k]][reverse[k]] ++;
                    }
                }

                if (m<2){
                    inter_dd[g][k] ++;
                    inter_dd[neigh[g][k]][reverse[k]] ++;

                    inter_dd[h][k] ++;
                    inter_dd[neigh[h][k]][reverse[k]] ++;
                }
            }

            if (m == 2){
                inter_dd[g][2] ++;
                inter_dd[neigh[g][2]][reverse[2]] ++;
                inter_dd[g][1] ++;
                inter_dd[neigh[g][1]][reverse[1]] ++;

                inter_dd[h][2] ++;
                inter_dd[neigh[h][2]][reverse[2]] ++;
                inter_dd[h][1] ++;
                inter_dd[neigh[h][1]][reverse[1]] ++;
            }

            if (m == 3){
                inter_dd[g][3] ++;
                inter_dd[neigh[g][3]][reverse[3]] ++;
                inter_dd[g][0] ++;
                inter_dd[neigh[g][0]][reverse[0]] ++;

                inter_dd[h][3] ++;
                inter_dd[neigh[h][3]][reverse[3]] ++;
                inter_dd[h][0] ++;
                inter_dd[neigh[h][0]][reverse[0]] ++;
            }
        }

        jc_new = 0;
        for (m=0;m<2*dim;m++){

            x = neigh[l][m];
            y = neigh[l_new][m];
            g = neigh_cd[l][m];
            h = neigh_cd[l_new][m];

            jc_new += ener_self(h);
            jc_new += ener_self(g);

            jc_new += ener_cd(l_new,m)/(double)(inter_cd[l_new][m]);

            if (cholpos[neigh_dc[g][m]] == 1){
                jc_new += ener_cd(neigh_dc[g][m],reverse[m])/(double)(inter_cd[neigh_dc[g][m]][reverse[m]]);
            }
            if (cholpos[neigh_dc[h][m]] == 1){
                jc_new += ener_cd(neigh_dc[h][m],reverse[m])/(double)(inter_cd[neigh_dc[h][m]][reverse[m]]);
            }

            for (k=0;k<2*dim;k++){

                if (cholpos[x] == 1){
                    jc_new += ener_cd(x,k)/(double)(inter_cd[x][k]);
                    if (cholpos[neigh[x][k]] == 1) jc_new += ener_cc(x,k)/(double)(inter_cc[x][k]);
                }

                if (cholpos[y] == 1){
                    jc_new += ener_cd(y,k)/(double)(inter_cd[y][k]);
                    if (cholpos[neigh[y][k]] == 1) jc_new += ener_cc(y,k)/(double)(inter_cc[y][k]);
                }

                if (m<2){
                    jc_new += ener_dd(g,k)/(double)(inter_dd[g][k]);
                    jc_new += ener_dd(h,k)/(double)(inter_dd[h][k]);
                }
            }

            if (m == 2){
                jc_new += ener_dd(g,2)/(double)(inter_dd[g][2]);
                jc_new += ener_dd(g,1)/(double)(inter_dd[g][1]);

                jc_new += ener_dd(h,2)/(double)(inter_dd[h][2]);
                jc_new += ener_dd(h,1)/(double)(inter_dd[h][1]);
            }

            if (m == 3){
                jc_new += ener_dd(g,3)/(double)(inter_dd[g][3]);
                jc_new += ener_dd(g,0)/(double)(inter_dd[g][0]);

                jc_new += ener_dd(h,3)/(double)(inter_dd[h][3]);
                jc_new += ener_dd(h,0)/(double)(inter_dd[h][0]);
            }
        }

        delta = jc_new - jc_old;

        if (exp(-delta/temperature) < random()*quotf){
            SWAP(cholpos[l],cholpos[l_new]);
        }
    }
}

// Counting-Functions

void count_neighbors()
{
    int i,j;
    int l,x;
    int m;

    for (i=1;i<=nr_types;i++){
        neightype_dc[i] = 0;
        for (j=1;j<=i;j++){
            neightype[i][j] = 0;
        }
    }
    neightype_cc = 0;

    for (l=0;l<nr_lattice;l++){
        for (m=0;m<2*dim;m+=2){
            x = neigh[l][m];
            if (particlet[l] >= particlet[x]){ 
                neightype[particlet[l]][particlet[x]] ++;
            }
            else{
                neightype[particlet[x]][particlet[l]] ++;
            }
        }
    }

    for (l=0;l<nr_lattice;l++){
        neightype_dc[particlet[l]] += count_dc(l);
        if (cholpos[l] == 1) neightype_cc += count_cc(i)/2;
    }
}

int count_dc(int l)
{
    int m,n=0;
    for (m=0;m<=3;m++) n+= cholpos[neigh_dc[l][m]];
    return (n);
}

int count_cc(int l)
{
    int m,n=0;
    for (m=0;m<=3;m++) n += cholpos[neigh[l][m]];
    return (n);
}

int count2_dd(int l,int m)
{
    int c,n=0;
    for (c=1;c<=cmax_dd;c++) n += cholpos[neigh2_dd[l][m][c]];
    return (n);
}

int count2_cd(int l,int m)
{
    int c,n=0;
    for (c=1;c<=cmax_cd;c++) n += cholpos[neigh2_cd[l][m][c]];
    return (n);
}

int count2_cc(int l,int m)
{
    int c,n=0;
    for (c=1;c<=cmax_cc;c++) n += cholpos[neigh2_cc[l][m][c]];
    return (n);
}

void count_order()
{
    int c,l;

    for (l=0;l<nr_lattice;l++){
        c = count_dc(l);
        p_MC[particlet[l]][order[l]] ++;
        p_MC_neigh[particlet[l]][order[l]][c] ++;
    }
}

void norm_order()
{
    int s,c,i;
    double p_MC_norm;
    double p_MC_norm_neigh[5];

    for (i=1;i<=nr_types;i++){
        p_MC_norm = 0;
        for (s=0;s<=150;s++){
            p_MC_norm += p_MC[i][s];
            for (c=0;c<=4;c++) p_MC_norm_neigh[c] += p_MC_neigh[1][s][c];
        }
        for (s=0;s<=150;s++){
            p_MC[i][s] /= (p_MC_norm);
            for (c=0;c<=4;c++){
                p_MC_neigh[i][s][c] /= p_MC_norm_neigh[c];
            }
        }
    }
}

// Calculating Energies

double calc_enertot()
{
    int l,m;
    double enertot;

    enertot = 0;
    for (l=0;l<nr_lattice;l++){
        for (m=0;m<2*dim;m++){
            enertot += ener_dd(l,m);
            if (cholpos[neigh_dc[l][m]] == 1){
                enertot += ener_cd(neigh_dc[l][m],reverse[m]);
            }
            if (cholpos[l] == 1){
                enertot += ener_cd(l,m);
                if (cholpos[neigh[l][m]] == 1){
                    enertot += ener_cc(l,m);
                }
            }
        }
    }
    enertot /= 2.0;
    for (l=0;l<nr_lattice;l++){
        enertot += ener_self(l);
    }
    return(enertot);
}

double ener_dd(int l, int m)
{
    int i = particlet[l];
    int x = neigh[l][m];
    int j = particlet[x];

    int c = count2_dd(l,m);
    int ci = count_dc(l);
    int cj = count_dc(x);


    return ((1/8.0)*(nneigh_dd[i][ci] + nneigh_dd[j][cj])*
            H_dd_avg[i][j][avg_order(l,x)][c]);
}

double ener_cd(int l, int m)
{
    int x = neigh_cd[l][m];
    int i = particlet[x];

    int c = count2_cd(l,m);
    int cm = count_cc(l);

    return ((1/4.0)*nneigh_cd[i][cm]*H_cd[i][order[x]][c]);
}

double ener_cc(int l, int m)
{
    int c = count2_cc(l,m);
    
    return (H_cc[c]);
}

double ener_self(int l)
{
    int i = particlet[l];
    double c = (double)count_dc(l);

    double ord = ((double)(order[l])-50.0)/100.0;
    return (H_self[i][order[l]] - temperature*((Entro[i][order[l]]+(b_S*ord*ord+a_S*ord)*(c+c_N*c*c))));
}

// Computational functions

void transfer1(int l, int* vec)
{
    int ll,m;

    ll = nr_lattice/length[1];
    for (m=1;m<=dim;m++){
        vec[m] = l/ll;
        l -= (ll*vec[m]);
        if (m<dim) ll /= length[m+1];
    }
}

int transfer2(int *vec)
{
    int m;
    int l = 0;
    int ll = 1;

    for (m=1;m<=dim;m++){
        l += vec[dim+1-m]*ll;
        ll *= length[dim+1-m];
    }
    return(l);
}

void sort(int n, double *arr, int *brr)
{
    int l,i;
    double a,b;

    for (l=2;l<=n;l++){
        a = arr[l];
        b = brr[l];
        i = l-1;
        while (i > 0 && arr[i] > a){
            arr[i+1] = arr[i];
            brr[i+1] = brr[i];
            i--;
        }
        arr[i+1] = a;
        brr[i+1] = b;
    }
}

void make_neighbors() 
{
  int i,j,k,l,vec[4],vecn[4],m,keep; 
  int c_1[4][10],c_2[4][10] ;

  nr_lattice = 1;
  for (i=1 ; i <= dim ; i++) 
  {
    nr_lattice *= length[i]; //number of lattice squares is product of lengths in all directions
  }
  nr_chol = (int)(c_chol * nr_lattice);
  if (nr_lattice > N) printf("Alarm \n");

  for (i=0 ; i < nr_lattice ; i++) 
  {
    transfer1(i,vec);
    l = 0;
    for (j=1 ; j<= dim ; j++)
    for (k=-1 ; k <= 1 ; k+=2)
    {
      keep = vec[j];
      vec[j] = per(vec[j] +  k,j);
      m = transfer2(vec);
      neigh[i][l] = m;   // counts the neighbors of site i, i.e. neigh[i][0],...,neigh[i][2*dim-1]
      vec[j] = keep;
      l++;
    }

    transfer1(i,vec);
    vecn[1] = per(vec[1] -1,1); vecn[2] = per(vec[2] -1,2); neigh_dc[i][0] = transfer2(vecn);
    vecn[1] = per(vec[1] ,1); vecn[2] = per(vec[2] ,2); neigh_dc[i][1] = transfer2(vecn);
    vecn[1] = per(vec[1] ,1); vecn[2] = per(vec[2] -1,2); neigh_dc[i][2] = transfer2(vecn);
    vecn[1] = per(vec[1] - 1 ,1); vecn[2] = per(vec[2] ,2); neigh_dc[i][3] = transfer2(vecn);

    transfer1(i,vec);
    vecn[1] = per(vec[1] ,1); vecn[2] = per(vec[2] ,2); neigh_cd[i][0] = transfer2(vecn);
    vecn[1] = per(vec[1]+1 ,1); vecn[2] = per(vec[2]+1 ,2); neigh_cd[i][1] = transfer2(vecn);
    vecn[1] = per(vec[1]+1 ,1); vecn[2] = per(vec[2] ,2); neigh_cd[i][2] = transfer2(vecn);
    vecn[1] = per(vec[1]  ,1); vecn[2] = per(vec[2]+1 ,2); neigh_cd[i][3] = transfer2(vecn);
     
  }
  // Starting point is always  (2,2); 

  cmax_dd = 6;
  c_1[0][1] = 0;  c_2[0][1] = 1;    
  c_1[0][2] = 0;  c_2[0][2] = 2;    
  c_1[0][3] = 1;  c_2[0][3] = 1;    
  c_1[0][4] = 1;  c_2[0][4] = 2;    
  c_1[0][5] = 2;  c_2[0][5] = 1;    
  c_1[0][6] = 2;  c_2[0][6] = 2;    

  c_1[1][1] = 1;  c_2[1][1] = 1;    
  c_1[1][2] = 1;  c_2[1][2] = 2;    
  c_1[1][3] = 2;  c_2[1][3] = 1;    
  c_1[1][4] = 2;  c_2[1][4] = 2;    
  c_1[1][5] = 3;  c_2[1][5] = 1;    
  c_1[1][6] = 3;  c_2[1][6] = 2;    

  c_1[2][1] = 1;  c_2[2][1] = 0;    
  c_1[2][2] = 2;  c_2[2][2] = 0;    
  c_1[2][3] = 1;  c_2[2][3] = 1;    
  c_1[2][4] = 2;  c_2[2][4] = 1;    
  c_1[2][5] = 1;  c_2[2][5] = 2;    
  c_1[2][6] = 2;  c_2[2][6] = 2;    

  c_1[3][1] = 1;  c_2[3][1] = 1;    
  c_1[3][2] = 2;  c_2[3][2] = 1;    
  c_1[3][3] = 1;  c_2[3][3] = 2;    
  c_1[3][4] = 2;  c_2[3][4] = 2;    
  c_1[3][5] = 1;  c_2[3][5] = 3;    
  c_1[3][6] = 2;  c_2[3][6] = 3;    


  for (i=0 ; i < nr_lattice ; i++) 
  { 
    for (j=0 ; j < 2*dim ; j++)
    {
      for (l = 1 ; l <= cmax_dd; l++)
      {
        transfer1(i,vec); 
        vecn[1] = per(c_1[j][l]+vec[1]-2,1);  
        vecn[2] = per(c_2[j][l]+vec[2]-2,2);  
        neigh2_dd[i][j][l] = transfer2(vecn);
      }
    }
  }

  cmax_cd = 5;
  c_1[0][1] = 1;  c_2[0][1] = 2;    
  c_1[0][2] = 2;  c_2[0][2] = 1;    
  c_1[0][3] = 3;  c_2[0][3] = 2;    
  c_1[0][4] = 2;  c_2[0][4] = 3;    
  c_1[0][5] = 1;  c_2[0][5] = 1;    

  c_1[1][1] = 1;  c_2[1][1] = 2;    
  c_1[1][2] = 2;  c_2[1][2] = 1;    
  c_1[1][3] = 3;  c_2[1][3] = 2;    
  c_1[1][4] = 2;  c_2[1][4] = 3;    
  c_1[1][5] = 3;  c_2[1][5] = 3;   

  c_1[2][1] = 1;  c_2[2][1] = 2;    
  c_1[2][2] = 2;  c_2[2][2] = 1;    
  c_1[2][3] = 3;  c_2[2][3] = 2;    
  c_1[2][4] = 2;  c_2[2][4] = 3;    
  c_1[2][5] = 3;  c_2[2][5] = 1;    

  c_1[3][1] = 1;  c_2[3][1] = 2;    
  c_1[3][2] = 2;  c_2[3][2] = 1;    
  c_1[3][3] = 3;  c_2[3][3] = 2;    
  c_1[3][4] = 2;  c_2[3][4] = 3;    
  c_1[3][5] = 1;  c_2[3][5] = 3;    
  for (i=0 ; i < nr_lattice ; i++) 
  { 
    for (j=0 ; j < 2*dim ; j++)
    {
      for (l = 1 ; l <= cmax_cd; l++)
      {
        transfer1(i,vec); 
        vecn[1] = per(c_1[j][l]+vec[1]-2,1);  
        vecn[2] = per(c_2[j][l]+vec[2]-2,2);  
        neigh2_cd[i][j][l] = transfer2(vecn);
      }
    }
  }


  cmax_cc = 6;
  c_1[0][1] = 0;  c_2[0][1] = 2;    
  c_1[0][2] = 1;  c_2[0][2] = 1;    
  c_1[0][3] = 1;  c_2[0][3] = 3;    
  c_1[0][4] = 2;  c_2[0][4] = 1;    
  c_1[0][5] = 2;  c_2[0][5] = 3;    
  c_1[0][6] = 3;  c_2[0][6] = 2;    

  c_1[1][1] = 1;  c_2[1][1] = 2;    
  c_1[1][2] = 2;  c_2[1][2] = 1;    
  c_1[1][3] = 2;  c_2[1][3] = 3;    
  c_1[1][4] = 3;  c_2[1][4] = 1;    
  c_1[1][5] = 3;  c_2[1][5] = 3;    
  c_1[1][6] = 4;  c_2[1][6] = 2;    

  c_1[2][1] = 2;  c_2[2][1] = 0;    
  c_1[2][2] = 1;  c_2[2][2] = 1;    
  c_1[2][3] = 3;  c_2[2][3] = 1;    
  c_1[2][4] = 1;  c_2[2][4] = 2;    
  c_1[2][5] = 3;  c_2[2][5] = 2;    
  c_1[2][6] = 2;  c_2[2][6] = 3;    

  c_1[3][1] = 2;  c_2[3][1] = 1;    
  c_1[3][2] = 1;  c_2[3][2] = 2;    
  c_1[3][3] = 3;  c_2[3][3] = 2;    
  c_1[3][4] = 1;  c_2[3][4] = 3;    
  c_1[3][5] = 3;  c_2[3][5] = 3;    
  c_1[3][6] = 2;  c_2[3][6] = 4;    


  for (i=0 ; i < nr_lattice ; i++) 
  { 
    for (j=0 ; j < 2*dim ; j++)
    {
      for (l = 1 ; l <= cmax_dd; l++)
      {
        transfer1(i,vec); 
        vecn[1] = per(c_1[j][l]+vec[1]-2,1);  
        vecn[2] = per(c_2[j][l]+vec[2]-2,2);  
        neigh2_cc[i][j][l] = transfer2(vecn);
      }
    }
  }
}

int per(int i, int d) 
{
   if (i >= length[d] ) i -= length[d]; 
   if (i < 0) i += length[d]; 
   return(i);
}

int min_order(int l, int o)
{
    if (order[l] < order[o]){
        return order[l];
    }
    else{
        return order[o];
    }
}

int max_order(int l, int o)
{
    if (order[l] > order[o]){
        return order[l];
    }
    else{
        return order[o];
    }
}

int avg_order(int l, int o)
{
    return (order[l]+order[o]);
}

double avg_part(int l)
{
    double lm = (double)l;

    return (1/(1+exp((-lm+120))));
}

void read_input()
{
    FILE *input;
    char line[100];

    input = fopen("input.txt", "r");
    if (input == NULL){
        printf("Error: 'input.txt' not found!\n");
    }

    while (fgets(line,sizeof(line),input))
    {
        if (strstr(line, "Dimension: ") != NULL){
            sscanf(line, "Dimension: %d", &dim);
        }
        if (strstr(line, "DLiPC: ") != NULL){
            sscanf(line, "DLiPC: %d", &dlipc);
        }
        if (strstr(line, "Move Type: ") != NULL){
            sscanf(line, "Move Type: %d", &mv_type);
        }
        if (strstr(line, "Project Name: ") != NULL){
            sscanf(line, "Project Name: %s", &project_name);
        }
        if (strstr(line, "Number of Lipids: ") != NULL){
            sscanf(line, "Number of Lipids: %d", &nr_types);
        }
        if (strstr(line, "Lattice Size: ") != NULL){
            sscanf(line, "Lattice Size: %d", &lattice_size);
        }
        if (strstr(line, "Number of Steps: ") != NULL){
            sscanf(line, "Number of Steps: %d", &t_total);
        }
        if (strstr(line, "Number of Runs: ") != NULL){
            sscanf(line, "Number of Runs: %d", &nr_average);
        }
        if (strstr(line, "Data Points: ") != NULL){
            sscanf(line, "Data Points: %d", &nr_count);
        }
        if (strstr(line, "Temperature: ") != NULL){
            sscanf(line, "Temperature: %lf", &read_temperature);
        }
        if (strstr(line, "Orderwidth: ") != NULL){
            sscanf(line, "Orderwidth: %d", &orderwidth);
        }
        if (strstr(line, "Concentration (CHOL): ") != NULL){
            sscanf(line, "Concentration (CHOL): %lf", &c_chol);
        }
        if (strstr(line, "Concentration (DLiPC): ") != NULL){
            sscanf(line, "Concentration (DLiPC): %lf", &c_lip[1]);
        }
        if (strstr(line, "Equilibrium: ") != NULL){
            sscanf(line, "Equilibrium: %lf", &pequil);
        }
    }
    fclose(input);
}

// Plot Results

void plot_info()
{
    int l,i,s;
    double enertot,contacts;
    double c_corr = c_chol/(1-c_chol);

    printf("\n");
    printf("Number of Chol: %d\n", nr_chol);
    printf("Average Order Parameter:");
    for (i=1;i<=nr_types;i++){
        double av = 0;
        double sum = 0;
        for (s=0;s<=150;s++){
            av += (s-50)/100.0*p_MC[i][s];
            sum += p_MC[i][s];
        }
        for (s=0;s<=150;s++){
            p_MC[i][s] /= sum;
        }
        printf(" %lf",av/sum);
    }
    printf("\n");

    printf("\n");
    enertot = calc_enertot();
    printf("Total Energie: %lf\n",enertot);
    
    printf("\n");
    contacts = 0;
    for (l=0;l<=nr_lattice;l++){
        if (cholpos[l] == 1){
            contacts += (double)count_cc(l);
        }
    }
    contacts /= pow(c_corr,2)*4*(double)(nr_lattice);
    printf("Av. CC-Contacts: %lf\n",contacts);
}

void plot_lattice(int n, char ufname[])
{
    FILE *fconfig;
    char fname[1024], str[100];
    int l,x,vec[4];
    double vec_chol[4];

    strcpy(fname, project_name);
    strcat(fname, "_lattice_");
    sprintf(str, "%d", n);
    strcat(fname, str);
    strcat(fname, ".dat");

    fconfig = fopen(fname, "w");
    for (vec[1]=0;vec[1]<length[1];vec[1]++){
        for (vec[2]=0;vec[2]<length[2];vec[2]++){
            x = transfer2(vec);
            fprintf(fconfig, "%d %d %d %d\n", vec[1],vec[2],particlet[x],order[x]);
        }
    }

    for (l=0;l<nr_lattice;l++){
        if (cholpos[l] == 1){
            transfer1(l,vec);
            vec_chol[1] = (double)(vec[1]) + 0.5;
            vec_chol[2] = (double)(vec[2]) + 0.5;
            if (vec_chol[1] > (double)(length[1])) vec_chol[1] -= (double)(length[1]);
            if (vec_chol[2] > (double)(length[2])) vec_chol[2] -= (double)(length[2]);
            fprintf(fconfig, "%lf %lf\n",vec_chol[1],vec_chol[2]);
        }
    }
    fclose(fconfig);
}

void plot_order(int n, char ufname[])
{
    FILE *fconfig;
    char fname[1024], str[100];
    int s,i;
    int c;

    if ((c_chol)/(1+c_chol) == 0.1) c = 10;
    if ((c_chol)/(1+c_chol) == 0.2) c = 20;
    if ((c_chol)/(1+c_chol) == 0.3) c = 30;


    strcpy(fname, project_name);
    strcat(fname, "_");
    sprintf(str, "%d", c);
    strcat(fname, str);
    strcat(fname, "_order_");
    sprintf(str, "%d", n);
    strcat(fname, str);
    strcat(fname, ".dat");

    fconfig = fopen(fname, "w");
    
    for (s=0;s<=150;s++){
        fprintf(fconfig,"%d ",s);
        for (i=1;i<=nr_types;i++){
            fprintf(fconfig,"%lf ",p_MC[i][s]);
        }
        fprintf(fconfig,"%lf \n",p_MD[s]);
    }
    fclose(fconfig);
}

void plot_order_neigh(int n, char ufname[])
{
    FILE *fconfig;
    char fname[1024], str[100];
    int s,c,i,x;

    if ((c_chol)/(1+c_chol) == 0.1) c = 10;
    if ((c_chol)/(1+c_chol) == 0.2) c = 20;
    if ((c_chol)/(1+c_chol) == 0.3) c = 30;

    strcpy(fname, project_name);
    strcat(fname, "_");
    sprintf(str, "%d", c);
    strcat(fname, str);
    strcat(fname, "_neigh_");
    sprintf(str, "%d", n);
    strcat(fname, str);
    strcat(fname, ".dat");

    fconfig = fopen(fname, "w");
    
    for (s=0;s<=150;s++){
        fprintf(fconfig,"%d ",s);
        for (i=1;i<=nr_types;i++){
            for (x=0;x<=4;x++){
                fprintf(fconfig,"%lf ",p_MC_neigh[i][s][x]);
            }
        }
        fprintf(fconfig,"\n");
    }
    fclose(fconfig);
}

void plot_Entro(int n, char ufname[])
{
    FILE *fconfig;
    char fname[1024], str[100];
    int s,c,i;

    strcpy(fname, project_name);
    strcat(fname, "_entro_");
    sprintf(str, "%d", n);
    strcat(fname, str);
    strcat(fname, ".dat");

    fconfig = fopen(fname, "w");
    
    for (s=0;s<=150;s++){
        fprintf(fconfig,"%d ",s);
        for (i=1;i<=nr_types;i++){
            fprintf(fconfig, "%lf %lf", Entro[1][s], init_origin[1][s]);
            }
        
        fprintf(fconfig,"\n");
    }
    fclose(fconfig);
}

// Fit-Parameters

double NN_DD_L(int neighbor)
{
   if (neighbor == 0)
    {
      return (5.103242510197753 - 0.005085225361896779*read_temperature);
    }
    if (neighbor == 1)
    {
      return (5.113452890355549 - 0.005357476571283622*read_temperature);
    }
    if (neighbor == 2)
    {
      return (5.041564491528647 - 0.005448578551134591*read_temperature);
    }
    if (neighbor == 3)
    {
      return (4.907827350161882 - 0.005455112431505342*read_temperature);
    }
    if (neighbor == 4)
    {
      return (4.437067046342299 - 0.004487121970234406*read_temperature);
    } 
}  

double NN_DD_P(int neighbor)
{
    if (neighbor == 0)
    {
      return (5.29864598362317 + (-1.046036898442689/ 
      (1 + exp(-0.25679510115746684 * (read_temperature - 324.0478762880865)))));
    }
    if (neighbor == 1)
    {
      return (5.291045158728617 + (-1.155530326754457 / 
      (1 + exp(-0.28882041884126536 * (read_temperature - 324.2842194308451)))));
    }
    if (neighbor == 2)
    {
      return (5.2916826656954115 + (-1.2845977425568815 / 
      (1 + exp(-0.28664779185612405 * (read_temperature - 323.9102671158312)))));
    }
    if (neighbor == 3)
    {
      return (5.241010005432148 + (-1.3888048522693892 / 
      (1 + exp(-0.2828101614730355 * (read_temperature - 323.715223040918)))));
    }
    if (neighbor == 4)
    {
      return (5.1776041414015666 + (-1.4834873657631906 / 
      (1 + exp(-0.27016244921668536 * (read_temperature - 323.4582319860251)))));
    }
}

double NN_CD_L(int neighbor)
{
    if (neighbor == 0)
    {
      return (6.380169972613336 - 0.006504015670675491*read_temperature);
    }
    if (neighbor == 1)
    {
      return (6.191307287280158 -0.006244334905148628*read_temperature);
    }
    if (neighbor == 2)
    {
      return (6.023217470971802 -0.006163443775914513*read_temperature);
    }
    if (neighbor == 3)
    {
      return (5.688301398563751 -0.005732889536537225*read_temperature);
    }
    if (neighbor == 4)
    {
      return (4.953838034411138 -0.0041862507153207894*read_temperature);
    }
}

double NN_CD_P(int neighbor)
{
    if (neighbor == 0)
    {
      return (6.149307668194549 + (-1.2123337139100219/ 
      (1 + exp(-0.279420884957051 * (read_temperature - 323.5622172239552)))));
    }
    if (neighbor == 1)
    {
      return (6.132031203205369 + (-1.3391102930661072 / 
      (1 + exp(-0.2848108023841632 * (read_temperature - 323.4614496375261)))));
    }
    if (neighbor == 2)
    {
      return (6.09118413190892 + (-1.431194943529606 / 
      (1 + exp(-0.2632957680847098 * (read_temperature - 323.22226054097797)))));
    }
    if (neighbor == 3)
    {
      return (6.0063999161659885 + (-1.5233708107172301 / 
      (1 + exp(-0.26322941684422435 * (read_temperature - 323.3131242555826)))));
    }
    if (neighbor == 4)
    {
      return (5.836996032158986 + (-1.4812685885334167	 / 
      (1 + exp(-0.23159639245144456 * (read_temperature - 322.96436681509624)))));
    }
}

double H_self_L(int ord)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    return (-198.2144241389774
    - 5.767114535779552*o 
    - 21.392325238025208*o*o
    + 126.763192813133*o*o*o
    - 168.20406177882367*o*o*o*o);
}

double H_self_P(int ord)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    return (-223.3148767446365
    + 11.182283360518415*o 
    - 24.886674151265726*o*o
    - 23.449813607231356*o*o*o
    + 9.052427099706513*o*o*o*o);
}

double H_DD_LL(int ord,int neighbor)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    if (neighbor == 0)
    {
      return ((-64.2385777507702)
      + (-85.0/(1+exp(-5.5152388615934305*(o-1.0)))));
    }
    if (neighbor == 1)
    {
      return ((-61.907988669194644)
      + (-85.0/(1+exp(-5.434468636915715*(o-1.0)))));
    }
    if (neighbor == 2)
    {
      return ((-59.51808320647846)
      + (-85.0/(1+exp(-5.7222222085142835*(o-1.0)))));
    }
    if (neighbor == 3)
    {
      return ((-56.80010184207919)
      + (-85.0/(1+exp(-6.046112173136968*(o-1.0)))));
    }
    if (neighbor == 4)
    {
      return ((-54.10623762430083)
      + (-85.0/(1+exp(-5.73206006922409*(o-1.0)))));
    }
    if (neighbor == 5)
    {
      return ((-52.25671633964124)
      + (-85.0/(1+exp(-5.7767043757695635*(o-1.0)))));
    }
    if (neighbor == 6)
    {
      return ((-49.57073396130829)
      + (-65.0/(1+exp(-5.917923875494118*(o-1.0)))));
    }
}

double H_DD_PP_avg(int ord,int neighbor)
{
    double p = (double)ord/2.0;
    double o = (p-50.0)/100.0;

    double m = o*o;

    if (o > (0.0)) return (-34.14372*m*o + 4.44322*m - 0.29837*o - 64.32286 + 3.0*neighbor);

    else return (-64.0 + 3.0*neighbor);
}

double H_DD_LL_avg(int ord,int neighbor)
{
    double p = (double)ord/2.0;
    double o = (p-50.0)/100.0;

    double m = o*o;

    if (o > (0.0)) return (-51.34502*m*o + 11.80470*m - 0.20345*o - 64.72701 + 2.5*neighbor);

    else return (-64.35 + 2.5*neighbor);
}

double H_DD_PP(int ord,int neighbor)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    if (neighbor == 0)
    {
      return ((-64.01505834779738)
      + (-85.0/(1+exp(-5.790501529508871*(o-1.0)))));
    }
    if (neighbor == 1)
    {
      return ((-62.140061908466556)
      + (-71.26511634433881/(1+exp(-5.738073754577456*(o-1.0)))));
    }
    if (neighbor == 2)
    {
      return ((-59.13452426429441)
      + (-67.70009843013932/(1+exp(-5.5794520397459575*(o-1.0)))));
    }
    if (neighbor == 3)
    {
      return ((-56.25342719515398)
      + (-67.07890880720473/(1+exp(-5.6396362558218724*(o-1.0)))));
    }
    if (neighbor == 4)
    {
      return ((-52.77151128708946)
      + (-70.34178084725163/(1+exp(-5.87616410257819*(o-1.0)))));
    }
    if (neighbor == 5)
    {
      return ((-49.71869099631439)
      + (-65.0/(1+exp(-6.040477808247997*(o-1.0)))));
    }
    if (neighbor == 6)
    {
      return ((-43.935239780557815)
      + (-76.95590795060127/(1+exp(-5.0*(o-1.0)))));
    }
}

double H_DD_LP(int ord,int neighbor)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    if (neighbor == 0)
    {
      return ((-64.00634092896922)
      + (-85.0/(1+exp(-5.3300771721508475*(o-1.0)))));
    }
    if (neighbor == 1)
    {
      return ((-61.50693262332394)
      + (-85.0/(1+exp(-5.8495526115388845*(o-1.0)))));
    }
    if (neighbor == 2)
    {
      return ((-59.08429503488477)
      + (-85.0/(1+exp(-6.234155020540408*(o-1.0)))));
    }
    if (neighbor == 3)
    {
      return ((-56.3906908602111)
      + (-68.39937312113554/(1+exp(-6.0663592559087*(o-1.0)))));
    }
    if (neighbor == 4)
    {
      return ((-52.790928365243936)
      + (-85.0/(1+exp(-5.627236498118878*(o-1.0)))));
    }
    if (neighbor == 5)
    {
      return ((-49.083301346156844)
      + (-85.0/(1+exp(-5.589864043019059*(o-1.0)))));
    }
    if (neighbor == 6)
    {
      return ((-47.150888458064365)
      + (-85.0/(1+exp(-5.621171915591369*(o-1.0)))));
    }
}

double H_CC_L(int neighbor)
{
    return (-22.63587131894127 + 1.4330068242577039*neighbor);
}

double H_CC_P(int neighbor)
{
    return (-21.406267414625454 + 1.4939079023403155*neighbor);
}

double H_CD_L(int ord,int neighbor)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    if (neighbor == 0)
    {
      double q = pow(((0.7539372160118417-1.5)/(o-1.5)),1.6592956112169561);
      return (36.483321544385326*q*(q-2.0));
    }
    if (neighbor == 1)
    {
      double q = pow(((0.7641800265641312-1.5)/(o-1.5)),1.6255609292682194);
      return (34.13691226888765*q*(q-2.0));
    }
    if (neighbor == 2)
    {
      double q = pow(((0.7699158189596492-1.5)/(o-1.5)),1.8291013645164085);
      return (33.68025196550303*q*(q-2.0));
    }
    if (neighbor == 3)
    {
      double q = pow(((0.7228854333900939-1.5)/(o-1.5)),2.4111347746661385);
      return (30.4277779459905*q*(q-2.0));
    }
    if (neighbor == 4)
    {
      double q = pow(((0.6740667747083389-1.5)/(o-1.5)),2.5);
      return (26.601605192448396*q*(q-2.0));
    }
    if (neighbor == 5)
    {
      double q = pow(((0.7419107199339535-1.5)/(o-1.5)),2.2050855814147825);
      return (26.949507831362077*q*(q-2.0));
    }
}

double H_CD_P(int ord,int neighbor)
{
    double p = (double)ord;
  double o = (p-50.0)/100;
  if (neighbor == 0)
  {
    double q = pow(((0.8015034812730596-1.5)/(o-1.5)),1.2104586776529458);
    return (34.162896996372645*q*(q-2.0));
  }
  if (neighbor == 1)
  {
    double q = pow(((0.7469858015529757-1.5)/(o-1.5)),1.4201657855084255);
    return (31.950030406520444*q*(q-2.0));
  }
  if (neighbor == 2)
  {
    double q = pow(((0.7529462715453211-1.5)/(o-1.5)),1.4532308249852621);
    return (29.894418246485863*q*(q-2.0));
  }
  if (neighbor == 3)
  {
    double q = pow(((0.7569774323440122-1.5)/(o-1.5)),1.3923986491593072);
    return (27.70235738706128*q*(q-2.0));
  }
  if (neighbor == 4)
  {
    double q = pow(((0.7613793423645616-1.5)/(o-1.5)),1.6463683453577225);
    return (25.6221899051367*q*(q-2.0));
  }
  if (neighbor == 5)
  {
    double q = pow(((0.8164971196090005-1.5)/(o-1.5)),1.3387806908751303);
    return (23.2981162468298*q*(q-2.0));
  }
}

double E_L(int ord)
{
    double entro[] = {-32.585663	,-30.866203	,-29.232518	,-27.681936	,-26.211767	,-24.819306	,-23.501839	,-22.256655	,-21.081049	,-19.972329	,-18.927823	,-17.944885	,-17.020896	,-16.153274	,-15.339477	,-14.577006	,-13.86341	,-13.196286	,-12.57329	,-11.992131	,-11.450581	,-10.946472	,-10.4777	,-10.042228	,-9.638087	,-9.263374	,-8.916259	,-8.59498	,-8.297849	,-8.023245	,-5.879377	,-5.774855	,-5.680293	,-5.587116	,-5.495373	,-5.40958	,-5.256765	,-5.124981	,-4.841459	,-4.879505	,-4.76921	,-4.672118	,-4.584211	,-4.367889	,-4.280518	,-4.185753	,-4.147191	,-4.068844	,-4.038433	,-3.994268	,-3.952651	,-3.897551	,-3.86577	,-3.823503	,-3.795968	,-3.759546	,-3.742374	,-3.714611	,-3.702783	,-3.690262	,-3.667644	,-3.66403	,-3.649278	,-3.645	,-3.644549	,-3.651427	,-3.655617	,-3.665521	,-3.682499	,-3.701393	,-3.724392	,-3.75123	,-3.783737	,-3.820218	,-3.864088	,-3.909537	,-3.960334	,-4.018678	,-4.082179	,-4.151228	,-4.2264	,-4.30541	,-4.394138	,-4.486408	,-4.58696	,-4.693939	,-4.806593	,-4.92647	,-5.054721	,-5.188704	,-5.332174	,-5.484733	,-5.641997	,-5.809017	,-5.987276	,-6.174251	,-6.366742	,-6.573127	,-6.789823	,-7.018619	,-7.259292	,-7.51496	,-7.782613	,-8.063208	,-8.362974	,-8.677154	,-9.008832	,-9.362476	,-9.73028	,-10.127916	,-10.547617	,-10.9955	,-11.468903	,-11.974876	,-12.512817	,-13.090008	,-13.629436	,-14.195804	,-14.805544	,-15.447538	,-16.133729	,-13.942953	,-14.676721	,-15.473254	,-16.337913	,-17.276442	,-18.294987	,-19.400117	,-20.598843	,-21.898641	,-23.307475	,-24.833817	,-26.486672	,-28.275603	,-30.210754	,-32.302876	,-34.563352	,-37.004226	,-39.638228	,-42.478806	,-45.54015	,-48.837225	,-52.3858	,-56.20248	,-60.304737	,-64.710944	,-69.440406	,-74.513399	,-79.951198	,-85.776118	,-92.011551};
    return entro[ord];
}

double E_P(int ord)
{
    double entro[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,-4.118937,-3.852361,-3.772671,-3.697232,-3.546854,-3.496151,-3.340507,-3.16266,-2.961434,-2.842948,-2.713178,-2.532163,-2.389425,-2.249913,-2.112968,-1.978765,-1.802305,-1.673868,-1.367041,-1.26247,-1.153913,-1.039294,-1.019863,-0.871165,-0.708574,-0.632807,-0.549904,-0.469855,-0.402968,-0.350183,-0.28764,-0.243385,-0.194656,-0.177376,-0.125448,-0.122756,-0.078265,-0.0648940000000002,-0.057086,-0.077941,0,-0.0793690000000002,-0.0952250000000001,-0.126744,-0.0963850000000002,-0.144892,-0.250412,-0.246676,-0.350902,-0.404134,-0.4623,-0.572485,-0.649316,-0.753291,-0.866856,-0.987921,-1.12086,-1.259962,-1.409906,-1.56693,-1.691124,-1.888587,-2.077237,-2.256544,-2.464971,-2.680907,-2.909033,-3.147398,-3.399947,-3.659804,-3.929231,-4.191412,-4.489291,-4.793285,-5.095529,-5.440075,-5.76721,-6.116733,-6.484325,-6.852125,-7.257355,-7.63348,-8.052163,-8.486288,-8.930542,-9.379204,-9.853483,-10.340207,-10.842307,-11.354913,-11.886182,-12.42949,-12.990063,-13.562656,-14.152452,-14.755464,-15.375649,-16.008366,-16.660715,-17.328951,-18.008317,-18.708756,-19.420442,-20.15205,-20.900657,-21.669096,-22.447564,-23.244881,-24.051944,-24.88605,-25.729152,-26.597938,-27.497828,-28.38938,-29.287956,-30.240541,-31.08282,-31.97013,-32.841195,-33.703116,-34.599741};
    return entro[ord];
}

// IBI

void MD_order()
{
    double norm_sum = 0.0;

    for (int s=0; s<=150; s++){
        double p = (double)s;
        double o = (p-50.0)/100.0;
        if (dlipc != 0){
            if (c_chol == 0.0){
                p_MD[s] = (exp(-0.14122 + 7.51277*o -9.36903*o*o -4.43679*o*o*o -97.86418*o*o*o*o +192.92704*o*o*o*o*o+19.37517*o*o*o*o*o*o -168.20577*o*o*o*o*o*o*o));
            }
            if (c_chol == 0.1){
                p_MD[s] = norm_dist(o,0.266242429883346,0.15584716148881);
            }
            if (c_chol == 0.2){
                p_MD[s] = norm_dist(o,0.303489659288194,0.155431284082244);
            }
            if (c_chol == 0.3){
                p_MD[s] = norm_dist(o,0.352622871221787,0.153619794800945);
            }
        }
        if (dlipc == 0){
            if (c_chol == 0.0){
                p_MD[s] = (exp(-0.9767356 + 8.69286553*o -12.7808724*o*o +12.12000201*o*o*o -21.41776641*o*o*o*o + 7.14478559*o*o*o*o*o));
            }
            if (c_chol == 0.1){
                p_MD[s] = norm_dist(o,0.404805087945088,0.183008055146023);
            }
            if (c_chol == 0.2){
                p_MD[s] = norm_dist(o,0.498731806873557,0.174224169351057);
            }
            if (c_chol == 0.3){
                p_MD[s] = norm_dist(o,0.591086345394284,0.15239698365142);
            }
        }
    }

    for (int s=0; s<=150; s++){
        norm_sum += p_MD[s];
    }

    for (int s=0; s<=150; s++){
      p_MD[s] /= (norm_sum);
    }
}

double norm_dist(double ord, double mean, double sigma){
    double vorfaktor = sqrt(2*3.14159265359)*sigma;
    double exponent = (-0.5)*(pow((ord-mean),2))*(1/pow(sigma,2));
    double e_funktion = exp(exponent);


    return e_funktion/vorfaktor;
}

double L_init(int ord)
{
    double temp = temperature*120;
    init_origin[1][ord] = (log(p_MD[ord])) + (1/(temp*0.00831448664))*(0.5*H_dd_avg[1][1][ord*2][0]*nneigh_dd[1][0] + H_self[1][ord]); 
    return init_origin[1][ord];
}

double P_init(int ord)
{
    double temp = temperature*120;
    init_origin[1][ord] = (log(p_MD[ord])) + (1/(temp*0.00831448664))*(0.5*H_dd_avg[1][1][ord*2][0]*nneigh_dd[1][0] + H_self[1][ord]); 
    return init_origin[1][ord];
    
}

void update_Entro(double alpha)
{
    double p_min, p_max;
    double omega_diff;
    for (int i=0; i<=150; i++){
        if (p_MC[1][i] != 0 /*&& i>=120*/){
          omega_diff = log((p_MD[i]/p_MC[1][i]))*alpha*2;
          Entro[1][i] += omega_diff;
        }
        if (p_MC[1][i] == 0 && p_MD[i] != 0){
          Entro[1][i] += omega_diff*alpha/10.0;
        }
    }

    norm_Entro();
}

void norm_Entro()
{
    double Entro_norm = 0;
    double Entro_max = Entro[1][0];
    for (int i=0; i<=150; i++){
        if (Entro_max < Entro[1][i]) Entro_max = Entro[1][i];
    }
    for (int i=0; i<=150;i++){
        Entro[1][i] -= Entro_max;
    }
    for (int i=0; i<=150; i++){
        Entro[1][i] = exp(Entro[1][i]);
    }
    for (int i=0; i<=150; i++){
        Entro_norm += Entro[1][i];
    }
    for (int i=0; i<=150; i++){
        Entro[1][i] /= Entro_norm;
    }
    for (int i=0; i<= 150; i++){
        Entro[1][i] = log(Entro[1][i]);
    }
}


double order_diff()
{
    double diff;
    for (int s=0; s<=150; s++){
      diff += (p_MD[s]-p_MC[1][s])*(p_MD[s]-p_MC[1][s]);
    }
    return diff;
}

double order_diff_neigh(int n)
{
    double diff;
    for (int s=0; s<=150; s++){
        diff += (p_MD[s]-p_MC_neigh[1][s][n])*(p_MD[s]-p_MC_neigh[1][s][n]);
    }
    return diff;
}

void fisher_yate_shuffle(int *array)
{
    int size = nr_lattice;
    for (int i = size - 1; i>0; i--){
        int j = rand() % (i+1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

double avg_order_mc()
{
    double av = 0;
    double sum = 0;

    for (int s=0;s<=150;s++){
        double q = (double)s;
        av += q*(double)p_MC[1][s];
        sum += (double)p_MC[1][s];
    }

    return (av/sum);
}

int max_order_mc()
{
    double max = 0.0;
    int n = 1;

    for (int s=0;s<=150;s++){
        if (p_MC[1][s] > max){
            max = p_MC[1][s];
            n = s;
        }
    }

    return n;
}

int max_order_md()
{
    double max = 0.0;
    int n = 1;

    for (int s=0;s<=150;s++){
        if (p_MD[s] > max){
            max = p_MD[s];
            n = s;
        }
    }

    return n;
}

double avg_order_neigh(int n)
{
    double av = 0;
    double sum = 0;

    for (int s=0;s<=150;s++){
        double q = (double)s;
        av += q*(double)p_MC_neigh[1][s][n];
        sum += (double)p_MC_neigh[1][s][n];
    }

    return (av/sum);
}

double avg_order_md()
{
    double av = 0;
    double sum = 0;

    for (int s=0;s<=150;s++){
        double q = (double)s;
        av += q*p_MD[s];
        sum += p_MD[s];
    }

    return (av/sum);
}