#include <stdio.h>
#include <iostream>


#include<string>
#include<fstream>
#include<time.h>

using std::cout;
using std::string;
using std::ifstream;

#define min(a, b) (((a) < (b)) ? (a) : (b))

//Parametros da GPU
#define threadspblock 1024

//Parametros da Entrada
string in_file_name = "Caso.txt";

//Parametros da Saida
string out_file_name = "Saida.csv";
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

__device__ long long int *d_Cn;//Combinacoes pré calculadas
void load_case (double* &E, int* &meas_plan, int* &UMs);

void show_completition_percentage(long long int n_analysed_combs);

//Funcs da GPU
__global__ void step1_enumeration(int* combs, long long int combs_first_id,int n_MU, int card,long long int n_combs_in_wave, int kmax);


int main()
{ 
    int* UMs;     //Unidades de Medicao [n_MU]
    int* meas_plan; //Plano de medição [n_meas x 7]
    double* E;     // Matriz Covariancia E [n_meas x n_meas]
    load_case(E,meas_plan,UMs);
    cout<<"nMUs: " << n_MU << "; nMeds: " <<n_meas<< "; kmax: " <<kmax<< '\n';

    
    // Conjunto Solucao
    //int n_sols=0;
    int* Sols;
    Sols = (int*)malloc(kmax * sol_Size * sizeof(int));

    //Combinacoes de elementos
    int* combs; 
    combs = (int*)malloc((size_t)wave_size * kmax * sizeof(int));

    // Vetor booleano  1: Combinacao critica 0: Combinacao nao
    int* is_crit; 
    is_crit = (int*)malloc(wave_size * sizeof(int));
    for(int i = 0; i<wave_size;i++) is_crit[i] = 1;

    //Alocacoes na GPU
	double *d_E;// Matrix de Covariancia
	cudaMalloc(&d_E,n_meas*n_meas * sizeof(double));
	cudaMemcpy(d_E,E,n_meas*n_meas * sizeof(double),cudaMemcpyHostToDevice);

    cudaMalloc(&d_Cn,(n_rows_Cn)*(n_colums_Cn)*sizeof( long long int));
    cudaMemcpyToSymbol(d_Cn,&Cn,(n_rows_Cn)*(n_colums_Cn)* sizeof(long long int),0,cudaMemcpyHostToDevice);
	int *d_combs;//Matriz com combinacoes enumeradas
    cudaMalloc(&d_combs,wave_size*kmax * sizeof(int));
	
    int *d_isCrit;//Matris que indica combinacoes criticas
	cudaMalloc(&d_isCrit,wave_size* sizeof(int));
	cudaMemset( d_isCrit,1,wave_size* sizeof(int));

    int* d_UMs;     //Unidades de Medicao [n_MU]
    cudaMalloc(&d_UMs,n_MU*sizeof(int));
    cudaMemcpy(d_UMs,UMs,n_MU*sizeof(int),cudaMemcpyHostToDevice);

    int* d_meas_plan; //Plano de medição [n_meas x 7]
    cudaMalloc(&d_meas_plan,n_meas*7*sizeof(int));
    cudaMemcpy(d_meas_plan,meas_plan,n_meas*7*sizeof(int),cudaMemcpyHostToDevice);

    // Numero de combinacoes avaliadas simultaneamentes
    int card;
    long long int n_combs_in_wave;

    for (card = 1; card <=kmax ;card++)
    {
        cout << "->Cardinalidade " << card<<'\n';
        cout << "Iniciado...\n";
        
        long long int n_analysed_combs = 0; // Combinacoes vizitadas em todas as ondas
        //FILE* combs_file;
        while (n_analysed_combs < Cn[n_MU * (n_colums_Cn) + card])
        {
            n_combs_in_wave = min(wave_size, Cn[n_MU * (n_colums_Cn) + card] - n_analysed_combs);
            //Enumerar
            step1_enumeration<<<wave_size/threadspblock,threadspblock>>>(d_combs,n_analysed_combs+1,n_MU,card,n_combs_in_wave,kmax);
            cudaMemcpy(combs,d_combs,wave_size*kmax * sizeof(int),cudaMemcpyDeviceToHost);
            for (int i = 0; i<n_combs_in_wave; i++)
            {
                for (int j = 0; j < kmax; j++)
                {
                    cout<<combs[i*kmax + j];
                }
                cout<<'\n';
            }
            //Prop. 2
            //Prop. 1
            //Atualização

            n_analysed_combs += wave_size;
            //Printa percentuais para acompanhar andamento da analise de criticalidades
            show_completition_percentage(n_analysed_combs);
            
        }
        
        printf("\nFinalizado!\n"); p25=0; p50=0;p75=0; 
        //fclose(combs_file);
    }
    
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

__global__ void step1_enumeration(int* combs, long long int combs_first_id, int n_MU, int card, long long int n_combs_in_wave, int kmax)
{
    int linha = threadIdx.x + blockDim.x*blockIdx.x;
    printf("%lli ", d_Cn[5]);
    // if(linha<n_combs_in_wave)
    // {  
    //     int nZ = n_MU - card;
    //     int nO = card;
    //     long long int n = linha + combs_first_id;
    //     for (long long int i = 0; i < n_MU; i++)
    //     {
    //         nZ--;
    //         if(nZ>=0){
    //             long long int zcomb = d_Cn[(n_MU - 1 - i) * (n_colums_Cn) +  min(nO,nZ)];
    //             if (zcomb < n)
    //             {
    //                 combs[linha * kmax + (card-nO)] = i;
    //                 nO--;
    //                 nZ++;
    //                 n = n - zcomb;
    //             }
    //         }else
    //         {
    //             combs[linha * kmax + (card-nO)] = i;
    //             nO --;
    //         }
    //     }
    //     for(int j=card; j<kmax;j++) {
    //          combs[linha * kmax + j] = -1;
    //     }
    // }
    
}  