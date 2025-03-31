#ifndef PARAMS_C
#define PARAMS_C

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

#include"../include/params.h"

void remove_white_lines_and_comments(FILE *input_fp)
{
	int temp_i;
	temp_i=getc(input_fp);
	if(temp_i=='\n' || temp_i==' ' || temp_i=='\043') // scan for white lines and comments
	{
		ungetc(temp_i, input_fp);

		temp_i=getc(input_fp);
		if(temp_i=='\n' || temp_i==' ') // white line
		{
			do
			{
				temp_i=getc(input_fp);
			}
			while(temp_i=='\n' || temp_i==' ');
		}
		ungetc(temp_i, input_fp);

		temp_i=getc(input_fp);
		if(temp_i=='\043')  // comment
		{
			do
			{
				temp_i=getc(input_fp);
			}
			while(temp_i!='\n');
		}
		else
		{
			ungetc(temp_i, input_fp);
		}

		remove_white_lines_and_comments(input_fp);
	}
	else
	{
		ungetc(temp_i, input_fp);
	}
}

void read_input(char *input_file_name, inv_params *param)
{
	FILE *input_fp;
	char str[STD_STRING_LENGTH], temp_str[STD_STRING_LENGTH];
	double temp_d;
	int temp_i, i;
    int err, end=1;
	input_fp=fopen(input_file_name, "r");  // open the input file
	if(input_fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{	
		while(end==1) // slide the file
		{
			remove_white_lines_and_comments(input_fp);
			err=fscanf(input_fp, "%s", str);
			if(err!=1)
			{
				fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			// lettura di Ns(S) e Nt(beta)
			if(strncmp(str, "size", 4)==0)
			{
				for(i=0; i<2; i++)
				{
					err=fscanf(input_fp, "%d", &temp_i);
					if(err!=1)
					{
						fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
						exit(EXIT_FAILURE);
					}
					if(i==0)
					{
						param->Ns=temp_i;
					}
					else
					{
						param->Nt=temp_i;
					}
				}
			}
			// lettura punti per l'inversione n_points(Nt)
			else if(strncmp(str, "points_for_inversion", 20)==0)
			{
				err=fscanf(input_fp, "%d", &temp_i);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            	param->n_points=temp_i;
			}
			// lettura punti scartati trash
			else if(strncmp(str, "trash", 5)==0)
			{
            	err=fscanf(input_fp, "%d", &temp_i);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            	param->trash=temp_i;
			}
			// lettura nome file con bootstrap
			else if(strncmp(str, "boostrap_corr_file", 18)==0)
			{
				err=fscanf(input_fp, "%s", temp_str);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
				strcpy(param->boot_corr_file, temp_str);
			}
			// lettura n_boot
			else if(strncmp(str, "n_boot", 6)==0)
			{
            	err=fscanf(input_fp, "%d", &temp_i);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            	param->Nboot=temp_i;
			}
			// lettura nome file con correlatore medio
			else if(strncmp(str, "mean_corr_file", 14)==0)
			{
				err=fscanf(input_fp, "%s", temp_str);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
				strcpy(param->mean_corr_file, temp_str);
			}
			// lettura Estar
        	else if(strncmp(str, "omega", 5)==0)
            { 
            	err=fscanf(input_fp, "%lf", &temp_d);
                if(err!=1)
            	{
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
            	}
            	param->Estar=temp_d;
            }
			// lettura apar
        	else if(strncmp(str, "apar", 4)==0)
            { 
            	err=fscanf(input_fp, "%lf", &temp_d);
                if(err!=1)
            	{
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
            	}
            	param->apar=temp_d;
            }
			// lettura Lside
        	else if(strncmp(str, "Lside", 5)==0)
            { 
            	err=fscanf(input_fp, "%lf", &temp_d);
                if(err!=1)
            	{
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
            	}
            	param->dLim=temp_d;
            }
			// lettura Rside
			else if(strncmp(str, "Rside", 5)==0)
            { 
            	err=fscanf(input_fp, "%lf", &temp_d);
                if(err!=1)
            	{
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
            	}
            	param->uLim=temp_d;
            }
			// lettura nome file output
			else if(strncmp(str, "output_file", 11)==0)
			{
				err=fscanf(input_fp, "%s", temp_str);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
				strcpy(param->out_file, temp_str);
			}
			// lettura nome file lambda vs rho
			else if(strncmp(str, "rho_vs_lambda_file", 18)==0)
			{
				err=fscanf(input_fp, "%s", temp_str);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
				strcpy(param->rho_vs_lambda_file, temp_str);
			}
			// lettura nome file delta_smear
			else if(strncmp(str, "smear_func_file", 15)==0)
			{
				err=fscanf(input_fp, "%s", temp_str);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
				strcpy(param->smear_delta_rec_file, temp_str);
			}
			// lettura nome file RS
			else if(strncmp(str, "rs_vs_lambda_file", 17)==0)
			{
				err=fscanf(input_fp, "%s", temp_str);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
				strcpy(param->rs_vs_lambda, temp_str);
			}
			// lettura sigma (leggo come stringa e poi converto in Real)
			else if(strncmp(str, "sigma", 5)==0)
			{
				err=fscanf(input_fp, "%s", temp_str);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
				Real temp_real_sigma(temp_str);
				param->sigma=temp_real_sigma;
			}
			// lettura Ezero (leggo come stringa e poi converto in Real)
			else if(strncmp(str, "Ezero", 5)==0)
			{
				err=fscanf(input_fp, "%s", temp_str);
                if(err!=1)
                {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
				Real temp_real_ezero(temp_str);
				param->E0=temp_real_ezero;				
			}
    		else
        	{
        		fprintf(stderr, "Error: unrecognized option %s in the file %s (%s, %d)\n", str, input_file_name, __FILE__, __LINE__);
        		exit(EXIT_FAILURE);
        	}

        	remove_white_lines_and_comments(input_fp);

       		// check if the read line is the last one
        	temp_i=getc(input_fp);
        	if(temp_i==EOF)
        	{
        		end=0;
        	}
        	else
        	{
            	ungetc(temp_i, input_fp);
        	}
        }
	}
    fclose(input_fp);
}

const int Nlambda = 180;  // Definition of Nlambda
PrecVec lambda = PrecVec::Zero(Nlambda);
Real lambda1, lambda2;
int ilambda1, ilambda2;


#if defined(HLN)
Real infLimit = 0.0;
#endif

#endif