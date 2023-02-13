#include <stdio.h>
#include <iostream>

#include<omp.h>
#include<math.h>

#include<string>
#include<fstream>
#include<time.h>

using namespace std;

#define min(a, b) (((a) < (b)) ? (a) : (b))

//Parametros da Entrada
string in_file_name = "Caso.txt";

//Parametros da Saida
string out_file_name = "Saida.csv";
string time_out_file_name = "Saida_tempos.csv";
int sol_Size = 10000;

//Parametros do Sitema
int n_MU; // Numero de Unidades de medicao no plano
int n_meas; // Numero de medidas do plano
int kmax; //Cardinalidade maxima avaliada

//Estutura da matriz de resultados de combinações (Cn)
string Cn_file_name = "combs1000em5.txt";
const int n_rows_Cn=1001;
const int n_colums_Cn=6;
long long int Cn[(n_rows_Cn) * (n_colums_Cn)] = { 0 };

// Numero de combinacoes avaliadas simultaneamentes
const long long int wave_size = (int)pow(2, 20);
int card;
long long int n_combs_in_wave;

int max_mat_size=0;

//Variaveis de percentual
bool p25 = 0;
bool p50 = 0;
bool p75 = 0;

int y=0;

//Variaveis de tempo
struct t_results{ 
    clock_t t_start, t_end;
    clock_t t_start_card[10], t_end_card[10];
};

void load_case (double* &E, int* &meas_plan, int* &UMs);

void step1_enumeration(int* combs, long long int combsIdini);

void step2_evaluation(int *is_crit, double* E, int* combs, int* MUs, int* meas_plan);
void build_aux_covariance(double *Ei, double* E, int* meas, int n_um_meas, int n_comb);
bool is_invertible(double *mat, int m);

void step3_confirmation(int* is_crit, int *combs, int* conjSol, int sol);

void step4_update(int* nSols, int* Sols, int* iscrit, int* combs);

void show_completition_percentage(long long int n_analysed_combs);

void save_results(int * Sols, int n_sols,t_results times, int* meas_plan);

int main()
{ 
    int* UMs;     //Unidades de Medicao [n_MU]
    int* meas_plan; //Plano de medição [n_meas x 7]
    double* E;     // Matriz Covariancia E [n_meas x n_meas]
    load_case(E,meas_plan,UMs);
    cout<<"nMUs: " << n_MU << "; nMeds: " <<n_meas<< "; kmax: " <<kmax<< '\n';


    t_results times;
    // Conjunto Solucao
    int nSols=0;
    int* Sols;
    Sols = (int*)malloc(kmax * sol_Size * sizeof(int));

    //Combinacoes de elementos
    int* combs; 
    combs = (int*)malloc((size_t)wave_size * kmax * sizeof(int));

    // Vetor booleano  1: Combinacao critica 0: Combinacao nao
    int* is_crit; 
    is_crit = (int*)malloc(wave_size * sizeof(int));
    times.t_start = clock();
    for (card = 1; card <=kmax ;card++)
    {
        cout << "->Cardinalidade " << card<<'\n';
        cout << "Iniciado...\n";
        times.t_start_card[card-1]=clock();
        long long int n_analysed_combs = 0; // Combinacoes vizitadas em todas as ondas
        //FILE* combs_file;
        while (n_analysed_combs < Cn[n_MU * (n_colums_Cn) + card])
        {
            n_combs_in_wave = min(wave_size, Cn[n_MU * (n_colums_Cn) + card] - n_analysed_combs);
            //Enumerar
            step1_enumeration(combs,n_analysed_combs+1);
            //Prop. 1
            step2_evaluation(is_crit, E, combs, UMs, meas_plan);
            int crits =0;
            //Prop. 2
            for (int sol =0; sol<nSols; sol++)
            {
                step3_confirmation(is_crit, combs, Sols, sol);
            }
            //Atualização
            step4_update(&nSols, Sols, is_crit, combs);

            n_analysed_combs += wave_size;
            
            //Printa percentuais para acompanhar andamento da analise de criticalidades
            show_completition_percentage(n_analysed_combs);
            
        }
        times.t_end_card[card-1] = clock();
        printf("\nFinalizado!\n"); p25=0; p50=0;p75=0; 
        //fclose(combs_file);
    }
    times.t_end = clock();

    save_results(Sols, nSols,times, UMs);

    free(UMs);
    free(meas_plan);
    free(E);
    free(combs);
    free(is_crit);
    free(Sols);
   
}

void load_case(double* &E, int* &meas_plan, int* &UMs)
{
    // Combinacoes pre calculadas (Cn)
    ifstream Cnfile(Cn_file_name.c_str());
	for (int i = 0; i < n_rows_Cn; i++)
	{
		for (int j = 0; j < n_colums_Cn; j++)
		{
			Cnfile >> Cn[i * (n_colums_Cn) + j];
		}
	}
	Cnfile.close();

    //Entrada programa
    ifstream in_file(in_file_name);
    
    in_file >> n_MU;
    in_file >> n_meas;
    in_file >> kmax;
    
    
    //Lista de unidades de medicao avaliadas
    UMs = (int*)malloc(n_MU*sizeof(int));
    for (int j = 0; j < n_MU; j++)
    {
        in_file >> UMs[j];
    }


    
    //Leitura do plano de med
    meas_plan = (int*)malloc(n_meas * 7 *sizeof(int));
    for (int i = 0; i < n_meas; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			in_file >> meas_plan[i * 7 + j];
		}
	}

    //Leitura da Matriz E
    E = (double*)malloc(n_meas * n_meas * sizeof(double));
    for (int i = 0; i < n_meas; i++)
	{
		for (int j = 0; j < n_meas; j++)
		{
			in_file >> E[i * n_meas + j];
		}
	}
}


void step1_enumeration(int* combs, long long int combsIdinicial)
{
    for (int linha= 0; linha< n_combs_in_wave; linha++){
        int nZ = n_MU - card;
        int nO = card;
        long long int n = linha + combsIdinicial;
        for (long long int i = 0; i < n_MU; i++)
        {
            nZ--;
            if(nZ>=0){
                long long int zcomb = Cn[(n_MU - 1 - i) * (n_colums_Cn) +  min(nO,nZ)];
                if (zcomb < n)
                {
                    combs[linha * kmax + (card-nO)] = i;
                    nO--;
                    nZ++;
                    n = n - zcomb;
                }
            }else
            {
                combs[linha * kmax + (card-nO)] = i;
                nO --;
            }
        }
        for(int j=card; j<kmax;j++) {
             combs[linha * kmax + j] = -1;
        }
    }
}  

                                                                                                                                                                                                                                                   
void step2_evaluation(int *is_crit, double* E, int* combs, int* MUs, int* meas_plan)
{
    
    for (int ind =0; ind<n_combs_in_wave; ind++){
        // Identifica medidas que fazem parte da UM
        int meas_number[50];
        int ind_meas = 0;
        for (int i = 0; i<card;i++)
        {
            int mu = MUs[combs[ind*kmax+i]];
            for (int meas =0; meas < n_meas; meas++)
            {
                if(meas_plan[meas*7 + 1] == mu) 
                {
                    meas_number[ind_meas] = meas_plan[meas*7];
                    ind_meas++;
                }
            }
        }
        if(max_mat_size<ind_meas) max_mat_size=ind_meas;
        // Avalia sub-matriz E
        double Ei[100*100];
        build_aux_covariance(Ei, E, meas_number, ind_meas, ind);
        
        if (!(is_invertible(Ei, ind_meas))) {
            is_crit[ind] = 1;
        }
        else {
            is_crit[ind] = 0;
		}
        free(meas_number);
        free(Ei);
    }
}

void build_aux_covariance(double *Ei, double* E, int* meas_number, int n_um_meas, int n_comb) 
{
    for (int i = 0; i < n_um_meas; i++)
	{
		for (int j = 0; j < n_um_meas; j++)
		{
			int m = meas_number[i]-1;
			int n = meas_number[j]-1;

			Ei[i * n_um_meas + j] = E[m * n_meas + n];
		}
	}
}

bool is_invertible(double *mat, int m) 
{
	bool inv;
	double pivo = 0.;
	for (int i = 0; i < m; i++) {

        //pivotiamento 
		int indmaior = i;
		double maior = mat[i * m + i];
		for (int j = i; j < m; j++)
		{
			if (abs(maior) < abs(mat[j * m + i]))
			{
				maior = mat[j * m + i];
				indmaior = j;
			}
		}
		for (int j = 0; j < m; j++)
		{
			double swap = mat[i * m + j];
			mat[i * m + j] = mat[indmaior * m + j];
			mat[indmaior * m + j] = swap;
		}
		pivo = mat[i*m + i];
		if (abs(pivo) < 0.0000000001) {
			inv = 0;
			return inv;
		}

		for (int j = 0; j < m; j++) {
			mat[i*m + j] = mat[i*m + j] / pivo;
		}


		for (int j = 0; j < m; j++) {
			if (j != i) 
            {
				pivo = mat[j*m + i];
				for (int l = 0; l < m; l++) {
					mat[j*m + l] = mat[j*m + l] - pivo * mat[i*m + l];
				}
			}
		}

	}

	inv = 1;
	return inv;
}


void step3_confirmation(int* is_crit, int *combs, int* conjSol, int sol)
{
    for (int crit = 0; crit < n_combs_in_wave; crit++)
    {
        if(is_crit[crit]==1)
        {
            int is = 0;
            int i = 0;
            int j = 0;
            while(conjSol[sol * kmax + i]!=-1 && i<card)
            {
                for (j = 0; j < card; j++)
                {
                    if (conjSol[sol * kmax + i]==combs[crit*kmax+j])
                        break;
                }
                if (j == card)
                    is = 1;
                i++;
            }
            is_crit[crit] = is;
        }
    }
}


void step4_update(int* nSols, int* Sols, int* iscrit, int* combs)
{
    for (int i = 0; i < n_combs_in_wave; i++)
    {
        if (iscrit[i] == 1)
        {
            for (int j = 0; j < kmax; j++)
                Sols[*nSols *kmax + j] = combs[i * kmax + j];
            *nSols += 1;
        }
    }
}


void show_completition_percentage(long long int n_analysed_combs)
{
    if (n_analysed_combs > Cn[n_MU * (n_colums_Cn) + card]/4 && p25==0)
    {
        printf("25%% ->");
        p25 = 1;
    }
    if (n_analysed_combs > Cn[n_MU * (n_colums_Cn) + card]/2 && p50==0)
    {
        printf("50%% ->");
        p50 = 1;
    }
    if (n_analysed_combs > Cn[n_MU * (n_colums_Cn) + card] * 3 / 4 && p75==0)
    {
        printf("75%%");
        p75 = 1;
    }
}

void save_results(int * Sols, int n_sols,t_results times, int* UMs)
{
    int n_crits[n_rows_Cn]={0};
    for (int i = 0; i<n_sols;i++){
        int j = 0; 
        while(Sols[i*kmax+j]!=-1 && j<kmax)
        {
            //cout << UMs[Sols[i*kmax+j]]<< ' ';
            j++;
        }
        n_crits[j-1]++;
        //cout << '\n';
    }

    double total_time = double(times.t_end - times.t_start) / double(CLOCKS_PER_SEC);
    printf( "tempo total: %f\n", total_time);

    double total_card_time[10];
    for ( int i = 0;  i<kmax; i++){
        total_card_time[i] = double(times.t_end_card[i] - times.t_start_card[i]) / double(CLOCKS_PER_SEC);   
        printf("Tempo total card %d: %f\n",i,total_card_time[i]);
    }

    //cout << "Maior maitriz invertida: "<< max_mat_size;
    FILE* Output_file;
    Output_file = fopen(out_file_name.c_str(), "w");
    fprintf(Output_file ,"Cardinalidade;numero de Cks\n");
    for (int i = 0; i < kmax; i++) {
        fprintf(Output_file ,"%i;%i\n", i + 1, n_crits[i]);
    }
    fprintf(Output_file ," Total de tuplas criticas de UM: %i\n", n_sols + 1);

    fprintf(Output_file, "Numero de combinacoes analisadas por cardinalidade\n");
    for (int i = 1; i <= kmax; i++)
        fprintf(Output_file, "%i; %lld\n", i, Cn[n_MU * (n_colums_Cn) + i]);
    
    
    fprintf(Output_file ,"Tempo total: \n");
    fprintf(Output_file ,"%f \n",total_time);
    fprintf(Output_file ,"Tempo total por cardinalidade: \n");
    for ( int i = 0;  i<kmax; i++){
        fprintf(Output_file ,"%d; %f\n",i+1,total_card_time[i]);
    }
   
    for (int i = 0; i<n_sols;i++){
        int j = 0; 
        while(Sols[i*kmax+j]!=-1 && j<kmax)
        {
            fprintf(Output_file,"UM%i;",  UMs[Sols[i*kmax+j]]);
            j++;
        }
        fprintf(Output_file,"\n");
    }
        
    
    fclose(Output_file);
}

