#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>

#define L 1 //1 2 3 4
#define LL 2 //2^(1 2 3 4)
#define LLL 4 //2^(2 3 4 5)

#define N 8

#define M 1 //# of trials
#define R 1 //the number of different initial conditions
#define S 100 //the number of different rules
#define T1 100 //transient
#define T2 1000 //time window used
int **mat;
int component[N],cnum;

void findcomponent(int node)
{
	int i;
	
	for (i=0;i<N;i++){
		if (mat[node][i]==1 && component[i]==0){
			component[i]=cnum;
			findcomponent(i);
		}
	}
}

double dabs(double x)
{
	if (x<0.0)return -x;
	else return x;
}

int main()
{
	FILE *fp,*fq,*fe,*fg,*fh,*fc,*fa;
	
	int i,j,k,l,t,m,r,s;
	int flag,flag1,flag2,tmp,tmp1,tmp2,tmp3,tmp4,ttmp;
	double dtmp;
	
	int x[N],xx[N][L+1],y[N];
	int in[N];
	int **nb,**rule;
	
	int check[N],sum,check_sum;
	
	double p1[N][LL],p2[N][LLL];
	double ****p3,****p4;
	int num1,num2,num3;
	
	double indeg[N],outdeg[N];
	
	double p=0.5;
	
	double q;
	int kav=2;
	int **adj,**adjflag;
	
	double *llap;
	double **lap,**edge,**grad,**circ;
	double div[N],sol[N];
	double diag[N];
	double **uu,**vv,**lap_inv;
	
	double el2_av,gl2_av,cl2_av,cratio_av;
	double el2,gl2,cl2,cratio;
	
	double *curl,*vp;
	int trinum;
	int **tri;
	
	double *ccurcur;
	double **curcur;
	double *diag2;
	double **uu2,**vv2,**curcur_inv;
	
	double **lcirc,**gcirc;
	
	double hl2_av,lcl2_av,gcratio_av,lcratio_av;
	double hl2,lcl2,gcratio,lcratio;
	
	double et,gt,ht,ct;
	double et_av,gt_av,ht_av,ct_av;
	
	char filename[100];
	
	gsl_matrix *v_gsl=gsl_matrix_alloc(N,N);
	gsl_vector *s_gsl=gsl_vector_alloc(N);
	gsl_vector *work=gsl_vector_alloc(N);
	gsl_matrix_view lap_gsl;
	
	gsl_matrix_view cur_gsl;
	
	
	const gsl_rng_type * TYPE;
	gsl_rng * ran;
	
	gsl_rng_env_setup();
	TYPE=gsl_rng_default;
	ran=gsl_rng_alloc(TYPE);
	gsl_rng_set(ran,2);
	
	
	adj=malloc(sizeof(int*)*N);
	adjflag=malloc(sizeof(int*)*N);
	
	mat=malloc(sizeof(int*)*N);
	lap=malloc(sizeof(double*)*N);
	edge=malloc(sizeof(double*)*N);
	grad=malloc(sizeof(double*)*N);
	circ=malloc(sizeof(double*)*N);
	lcirc=malloc(sizeof(double*)*N);
	gcirc=malloc(sizeof(double*)*N);
	
	
	uu=malloc(sizeof(double*)*N);
	vv=malloc(sizeof(double*)*N);
	lap_inv=malloc(sizeof(double*)*N);
	
	llap=malloc(sizeof(double)*(N*N));
	
	p3=malloc(sizeof(double***)*N);
	p4=malloc(sizeof(double***)*N);
	
	for (i=0;i<N;i++){
		adj[i]=malloc(sizeof(int)*N);
		adjflag[i]=malloc(sizeof(int)*N);
		
		mat[i]=malloc(sizeof(int)*N);
		
		lap[i]=malloc(sizeof(double)*N);
		edge[i]=malloc(sizeof(double)*N);
		grad[i]=malloc(sizeof(double)*N);
		circ[i]=malloc(sizeof(double)*N);
		lcirc[i]=malloc(sizeof(double)*N);
		gcirc[i]=malloc(sizeof(double)*N);
		
		uu[i]=malloc(sizeof(double)*N);
		vv[i]=malloc(sizeof(double)*N);
		lap_inv[i]=malloc(sizeof(double)*N);
		
		p3[i]=malloc(sizeof(double**)*N);
		p4[i]=malloc(sizeof(double**)*N);
		for (j=0;j<N;j++){
			p3[i][j]=malloc(sizeof(double*)*LL);
			p4[i][j]=malloc(sizeof(double*)*LLL);
			for (k=0;k<LL;k++){
				p3[i][j][k]=malloc(sizeof(double)*LL);
			}
			for (k=0;k<LLL;k++){
				p4[i][j][k]=malloc(sizeof(double)*LL);
			}
		}
	}
	
	sprintf(filename,"rtn_hodge_sw_av_n%d_k%d.txt",N,kav);
	fp=fopen(filename,"w");
	
	sprintf(filename,"rtn_hodge_sw_trial_n%d_k%d.txt",N,kav);
	fq=fopen(filename,"w");
	
	for (q=0.1;q<0.11;q=q+0.1){
		
		el2_av=0.0; gl2_av=0.0; cl2_av=0.0;
		cratio_av=0.0;
		
		hl2_av=0.0; lcl2_av=0.0;
		gcratio_av=0.0; lcratio_av=0.0;
		
		et_av=0.0; gt_av=0.0; ht_av=0.0; ct_av=0.0;
		
		for (m=0;m<M;m++){
			
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					adj[i][j]=0;
					adjflag[i][j]=0;
				}
			}
			for (i=0;i<N;i++){
				for (k=-kav;k<kav+1;k++){
					if (k!=0){
						adj[(i+k+N)%N][i]=1;
						adj[i][(i+k+N)%N]=1;
					}
				}
			}
			
			for (i=0;i<N;i++){
				for (k=-kav;k<kav+1;k++){
					if (k!=0 && gsl_rng_uniform(ran)<q){
						adjflag[(i+k+N)%N][i]=1;
						adjflag[i][(i+k+N)%N]=1;
					}
				}
			}
			
			for (i=0;i<N;i++){
				for (j=i+1;j<N;j++){
					if (adjflag[i][j]==1){
						if (gsl_rng_uniform(ran)<0.5){
							while (1){
								tmp=(int)((double)N*gsl_rng_uniform(ran));
								
								if (j!=tmp && adj[tmp][j]==0){
									adj[i][j]=0;
									adj[j][i]=0;
									adj[tmp][j]=1;
									adj[j][tmp]=1;
									break;
								}
							}
						} else{
							while (1){
								tmp=(int)((double)N*gsl_rng_uniform(ran));
								
								if (i!=tmp && adj[i][tmp]==0){
									adj[i][j]=0;
									adj[j][i]=0;
									adj[tmp][i]=1;
									adj[i][tmp]=1;
									break;
								}
							}
						}
					}
				}
			}
			
			for (i=0;i<N;i++){
				for (j=i+1;j<N;j++){
					if (adj[i][j]==1){
						if (gsl_rng_uniform(ran)<0.5)adj[i][j]=0;
						else adj[j][i]=0;
					}
				}
			}
			
			for (i=0;i<N;i++){
				in[i]=0;
				for (j=0;j<N;j++){
					in[i]+=adj[j][i];
				}
			}
			
			nb=malloc(sizeof(int*)*N);
			rule=malloc(sizeof(int*)*N);
			for (i=0;i<N;i++){
				nb[i]=malloc(sizeof(int)*in[i]);
				rule[i]=malloc(sizeof(int)*in[i]);
			}
			
			//in-neighbors
			for (i=0;i<N;i++){
				k=0;
				for (j=0;j<N;j++){
					if (adj[j][i]==1){
						nb[i][k]=j;
						k++;
					}
				}
			}
			//
			
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					edge[i][j]=0.0;
				}
			}
			
			for (s=0;s<S;s++){
				//random generation of rules
				for (i=0;i<N;i++){
					for (j=0;j<in[i];j++){
					//	rule[i][j]=gsl_ran_gaussian(ran,sigma);
						if (gsl_rng_uniform(ran)<0.5)rule[i][j]=1;
						else rule[i][j]=-1;
					}
				}
				
				for (k=0;k<N;k++){
					for (i=0;i<LL;i++){
						p1[k][i]=0.0;
					}
					
					for (i=0;i<LLL;i++){
						p2[k][i]=0.0;
					}
					
					for (l=0;l<in[k];l++){
						for (i=0;i<LL;i++){
							for (j=0;j<LL;j++){
								p3[k][nb[k][l]][i][j]=0.0;
							}
						}
						for (i=0;i<LLL;i++){
							for (j=0;j<LL;j++){
								p4[k][nb[k][l]][i][j]=0.0;
							}
						}
					}
				}
				
				for (r=0;r<R;r++){
					
					//initial state
					for (i=0;i<N;i++){
						if (gsl_rng_uniform(ran)<p)x[i]=1;
						else x[i]=-1;
					}
					
					//transient
					for (t=0;t<T1;t++){
						
						for (i=0;i<N;i++){
							tmp=0;
							for (j=0;j<in[i];j++){
								tmp+=rule[i][j]*x[nb[i][j]];
							}
							if (tmp<0)y[i]=-1;
							else y[i]=1;
						}
						
						for (i=0;i<N;i++){
							x[i]=y[i];
						}
					}
					
					for (i=0;i<N;i++){
						for (j=0;j<L;j++){
							xx[i][j]=0;
						}
						xx[i][L]=x[i];
					}
					
					//time window used to calculation
					for (t=0;t<T2+L-1;t++){
						for (i=0;i<N;i++){
							tmp=0;
							for (j=0;j<in[i];j++){
								tmp+=rule[i][j]*x[nb[i][j]];
							}
							if (tmp<0)y[i]=-1;
							else y[i]=1;
							
							for (j=0;j<L;j++){
								xx[i][j]=xx[i][j+1];
							}
							xx[i][L]=y[i];
						}
						
						if (t>=L-1){
							for (i=0;i<N;i++){
								num1=0; num2=0;
								for (k=0;k<L;k++){
									num1+=((xx[i][k]+1)/2)*(int)pow(2.0,k);
									num2+=((xx[i][k]+1)/2)*(int)pow(2.0,k);
								}
								num2+=((xx[i][L]+1)/2)*(int)pow(2.0,L);
								
								p1[i][num1]+=1.0;
								p2[i][num2]+=1.0;
								
								for (j=0;j<in[i];j++){
									num3=0;
									for (k=0;k<L;k++){
										num3+=((xx[nb[i][j]][k]+1)/2)*(int)pow(2.0,k);
									}
									
									p3[i][nb[i][j]][num1][num3]+=1.0;
									p4[i][nb[i][j]][num2][num3]+=1.0;
								}
							}
						}
						
						for (i=0;i<N;i++){
							x[i]=y[i];
						}
					}
				}//end of r-loop
				
				for (k=0;k<N;k++){
					for (i=0;i<LL;i++){
						p1[k][i]=p1[k][i]/(double)(T2*R);
					}
					
					for (i=0;i<LLL;i++){
						p2[k][i]=p2[k][i]/(double)(T2*R);
					}
					
					for (l=0;l<in[k];l++){
						for (i=0;i<LL;i++){
							for (j=0;j<LL;j++){
								p3[k][nb[k][l]][i][j]=p3[k][nb[k][l]][i][j]/(double)(T2*R);
							}
						}
						for (i=0;i<LLL;i++){
							for (j=0;j<LL;j++){
								p4[k][nb[k][l]][i][j]=p4[k][nb[k][l]][i][j]/(double)(T2*R);
							}
						}
					}
				}
				
				//construction of edge flow
				for (i=0;i<N;i++){
					for (j=0;j<in[i];j++){
						for (num2=0;num2<LLL;num2++){
							for (num3=0;num3<LL;num3++){
								if (p4[i][nb[i][j]][num2][num3]>0.0)edge[nb[i][j]][i]+=p4[i][nb[i][j]][num2][num3]*(log(p4[i][nb[i][j]][num2][num3])/log(2.0) + log(p1[i][num2%LL])/log(2.0) - log(p2[i][num2])/log(2.0) - log(p3[i][nb[i][j]][num2%LL][num3])/log(2.0));
							}
						}
					}
				}
			}//end of s-loop
			
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					edge[i][j]=edge[i][j]/(double)S;
				}
			}
			
			for (i=0;i<N;i++){
				for (j=i+1;j<N;j++){
					dtmp=edge[i][j]-edge[j][i];
					edge[i][j]=dtmp;
					edge[j][i]=-dtmp;
				}
			}
			
			//finding the number of components
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					mat[i][j]=0;
				}
			}
			
			for (i=0;i<N;i++){
				for (j=0;j<in[i];j++){
					if (i!=nb[i][j]){
						mat[nb[i][j]][i]=1;
						mat[i][nb[i][j]]=1;
					}
				}
			}
			
			tmp=0;
			for (i=0;i<N;i++){
				for (j=i+1;j<N;j++){
					if (mat[i][j]==1)tmp++;
				}
			}
			et=(double)tmp;
			
			for (i=0;i<N;i++){
				component[i]=0;
			}
			
			cnum=1; component[0]=cnum;
			findcomponent(0);
			for (i=1;i<N;i++){
				if (component[i]==0){
					cnum++;
					component[i]=cnum;
					findcomponent(i);
				}
			}
			
			gt=(double)(N-cnum);
			
			//graph Laplacian
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					lap[i][j]=0.0;
				}
			}
			
			for (i=0;i<N;i++){
				lap[i][i]=0.0;
				for (j=0;j<N;j++){
					lap[i][i]+=(double)mat[i][j];
				}
			}
			
			for (i=0;i<N;i++){
				for (j=0;j<in[i];j++){
					if (i!=nb[i][j]){
						lap[nb[i][j]][i]=-1.0;
						lap[i][nb[i][j]]=-1.0;
					}
				}
			}
			
			//generalized inverse of lap
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					llap[i*N+j]=lap[i][j];
				}
			}
			
			lap_gsl=gsl_matrix_view_array(llap,N,N);
			gsl_linalg_SV_decomp(&lap_gsl.matrix,v_gsl,s_gsl,work);
			
			for (i=0;i<N-cnum;i++){
				diag[i]=1.0/gsl_vector_get(s_gsl,i);
			}
			for (i=N-cnum;i<N;i++){
				diag[i]=0.0;
			}
			
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					uu[i][j]=gsl_matrix_get(&lap_gsl.matrix,i,j);
					vv[i][j]=gsl_matrix_get(v_gsl,i,j);
				}
			}
			
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					lap_inv[i][j]=0.0;
					for (k=0;k<N;k++){
						lap_inv[i][j]+=vv[i][k]*diag[k]*uu[j][k];
					}
				}
			}
			//
			
			
			//divergence of edge flow
			for (i=0;i<N;i++){
				div[i]=0.0;
				for (j=0;j<N;j++){
					div[i]+=edge[i][j];
				}
			}
			
			//solution (negative potential)
			for (i=0;i<N;i++){
				sol[i]=0.0;
				for (j=0;j<N;j++){
					sol[i]-=lap_inv[i][j]*div[j];
				}
			}
			
			//gradient flow
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					grad[i][j]=sol[j]-sol[i];
				}
			}
			
			//circular flow
			for (i=0;i<N;i++){
				for (j=0;j<N;j++){
					circ[i][j]=edge[i][j]-grad[i][j];
				}
			}
			
			el2=0.0; gl2=0.0; cl2=0.0;
			for (i=0;i<N;i++){
				for (j=i+1;j<N;j++){
					if (lap[i][j]<0.0){
						el2+=edge[i][j]*edge[i][j];
						gl2+=grad[i][j]*grad[i][j];
						cl2+=circ[i][j]*circ[i][j];
					}
				}
			}
			
			if (el2>0.0)cratio=cl2/el2;
			else cratio=0.0;
			
			el2_av+=el2/(double)M;
			gl2_av+=gl2/(double)M;
			cl2_av+=cl2/(double)M;
			cratio_av+=cratio/(double)M;
			
			//decomposing circular flow into harmonic and curl flows
			trinum=0;
			for (i=0;i<N;i++){
				for (j=i+1;j<N;j++){
					for (k=j+1;k<N;k++){
						if (mat[i][j]*mat[j][k]*mat[k][i]!=0)trinum++;
					}
				}
			}
			
			if (trinum>0){
				curl=malloc(sizeof(double)*trinum);
				tri=malloc(sizeof(int*)*trinum);
				for (i=0;i<trinum;i++){
					tri[i]=malloc(sizeof(int)*3);
				}
				
				l=0;
				for (i=0;i<N;i++){
					for (j=i+1;j<N;j++){
						for (k=j+1;k<N;k++){
							if (mat[i][j]*mat[j][k]*mat[k][i]!=0){
								tri[l][0]=i; tri[l][1]=j; tri[l][2]=k;
								curl[l]=circ[i][j]+circ[j][k]-circ[i][k];
								l++;
							}
						}
					}
				}
				
				ccurcur=malloc(sizeof(double)*trinum*trinum);
				diag2=malloc(sizeof(double)*trinum);
				vp=malloc(sizeof(double)*trinum);
				
				curcur=malloc(sizeof(double*)*trinum);
				uu2=malloc(sizeof(double*)*trinum);
				vv2=malloc(sizeof(double*)*trinum);
				curcur_inv=malloc(sizeof(double*)*trinum);
				for (i=0;i<trinum;i++){
					curcur[i]=malloc(sizeof(double)*trinum);
					uu2[i]=malloc(sizeof(double)*trinum);
					vv2[i]=malloc(sizeof(double)*trinum);
					curcur_inv[i]=malloc(sizeof(double)*trinum);
				}
				
				for (i=0;i<trinum;i++){
					for (j=0;j<trinum;j++){
						curcur[i][j]=0.0;
						if (tri[i][0]==tri[j][0] && tri[i][1]==tri[j][1])curcur[i][j]+=1.0;
						if (tri[i][0]==tri[j][0] && tri[i][1]==tri[j][2])curcur[i][j]-=1.0;
						if (tri[i][0]==tri[j][1] && tri[i][1]==tri[j][2])curcur[i][j]+=1.0;
						if (tri[i][1]==tri[j][0] && tri[i][2]==tri[j][1])curcur[i][j]+=1.0;
						if (tri[i][1]==tri[j][0] && tri[i][2]==tri[j][2])curcur[i][j]-=1.0;
						if (tri[i][1]==tri[j][1] && tri[i][2]==tri[j][2])curcur[i][j]+=1.0;
						if (tri[i][0]==tri[j][0] && tri[i][2]==tri[j][1])curcur[i][j]-=1.0;
						if (tri[i][0]==tri[j][0] && tri[i][2]==tri[j][2])curcur[i][j]+=1.0;
						if (tri[i][0]==tri[j][1] && tri[i][2]==tri[j][2])curcur[i][j]-=1.0;
					}
				}
				
				gsl_matrix *v_gsl2=gsl_matrix_alloc(trinum,trinum);
				gsl_vector *s_gsl2=gsl_vector_alloc(trinum);
				gsl_vector *work2=gsl_vector_alloc(trinum);
				
				//generalized inverse of curcur
				for (i=0;i<trinum;i++){
					for (j=0;j<trinum;j++){
						ccurcur[i*trinum+j]=curcur[i][j];
					}
				}
				
				cur_gsl=gsl_matrix_view_array(ccurcur,trinum,trinum);
				gsl_linalg_SV_decomp(&cur_gsl.matrix,v_gsl2,s_gsl2,work2);
				
				tmp=0;
				for (i=0;i<trinum;i++){
					dtmp=gsl_vector_get(s_gsl2,i);
					if (dtmp>1e-8)diag2[i]=1.0/dtmp;
					else {
						diag2[i]=0.0;
						tmp++;
					}
				}
				
				ct=(double)(trinum-tmp);
				ht=et-gt-ct;
				
				for (i=0;i<trinum;i++){
					for (j=0;j<trinum;j++){
						uu2[i][j]=gsl_matrix_get(&cur_gsl.matrix,i,j);
						vv2[i][j]=gsl_matrix_get(v_gsl2,i,j);
					}
				}
				
				for (i=0;i<trinum;i++){
					for (j=0;j<trinum;j++){
						curcur_inv[i][j]=0.0;
						for (k=0;k<trinum;k++){
							curcur_inv[i][j]+=vv2[i][k]*diag2[k]*uu2[j][k];
						}
					}
				}
				
				for (i=0;i<trinum;i++){
					vp[i]=0.0;
					for (j=0;j<trinum;j++){
						vp[i]+=curcur_inv[i][j]*curl[j];
					}
				}
				
				for (i=0;i<N;i++){
					lcirc[i][i]=0.0;
					for (j=i+1;j<N;j++){
						lcirc[i][j]=0.0;
						if (mat[i][j]!=0){
							for (k=0;k<trinum;k++){
								if (tri[k][0]==i && tri[k][1]==j)lcirc[i][j]+=vp[k];
								if (tri[k][0]==i && tri[k][2]==j)lcirc[i][j]-=vp[k];
								if (tri[k][1]==i && tri[k][2]==j)lcirc[i][j]+=vp[k];
							}
						}
					}
				}
				
				for (i=0;i<N;i++){
					for (j=0;j<i;j++){
						lcirc[i][j]=-lcirc[j][i];
					}
				}
				
				for (i=0;i<N;i++){
					for (j=0;j<N;j++){
						gcirc[i][j]=circ[i][j]-lcirc[i][j];
					}
				}
				
				hl2=0.0; lcl2=0.0;
				for (i=0;i<N;i++){
					for (j=i+1;j<N;j++){
						if (lap[i][j]<0.0){
							hl2+=gcirc[i][j]*gcirc[i][j];
							lcl2+=lcirc[i][j]*lcirc[i][j];
						}
					}
				}
				
				fe=fopen("edgeflow.txt","w");
				fg=fopen("gradientflow.txt","w");
				fh=fopen("harmonicflow.txt","w");
				fc=fopen("curlflow.txt","w");
				fa=fopen("adjmatrix.txt","w");
				for (i=0;i<N;i++){
					for (j=0;j<N;j++){
						if (adj[i][j]==1)fprintf(fa,"1 ");
						else fprintf(fa,"0 ");
						
						if (lap[i][j]<0.0){
							fprintf(fe,"%lf ",edge[i][j]);
							fprintf(fg,"%lf ",grad[i][j]);
							fprintf(fh,"%lf ",gcirc[i][j]);
							fprintf(fc,"%lf ",lcirc[i][j]);
						} else{
							fprintf(fe,"%lf ",0.0);
							fprintf(fg,"%lf ",0.0);
							fprintf(fh,"%lf ",0.0);
							fprintf(fc,"%lf ",0.0);
						}
					}
					fprintf(fe,"\n");
					fprintf(fg,"\n");
					fprintf(fh,"\n");
					fprintf(fc,"\n");
					fprintf(fa,"\n");
				}
				fclose(fe);
				fclose(fg);
				fclose(fh);
				fclose(fc);
				fclose(fa);
				
				
				gsl_matrix_free(v_gsl2);
				gsl_vector_free(s_gsl2);
				gsl_vector_free(work2);
				
				free(ccurcur);
				free(diag2);
				free(vp);
				
				if (curcur){
					for (i=0;i<trinum;i++){
						free(curcur[i]);
					}
				}
				if (uu2){
					for (i=0;i<trinum;i++){
						free(uu2[i]);
					}
				}
				if (vv2){
					for (i=0;i<trinum;i++){
						free(vv2[i]);
					}
				}
				if (curcur_inv){
					for (i=0;i<trinum;i++){
						free(curcur_inv[i]);
					}
				}
				
			} else{
				ct=0.0;
				ht=et-gt-ct;
				
				for (i=0;i<N;i++){
					for (j=0;j<N;j++){
						gcirc[i][j]=circ[i][j];
					}
				}
				
				hl2=0.0; lcl2=0.0;
				for (i=0;i<N;i++){
					for (j=i+1;j<N;j++){
						if (lap[i][j]<0.0){
							hl2+=gcirc[i][j]*gcirc[i][j];
						}
					}
				}
			}
			
			if (cl2>0.0)gcratio=hl2/el2;
			else gcratio=0.0;
			
			hl2_av+=hl2/(double)M;
			lcl2_av+=lcl2/(double)M;
			gcratio_av+=gcratio/(double)M;
			
			lcratio=cratio-gcratio;
			lcratio_av+=lcratio/(double)M;
			
			et_av+=et/(double)M;
			gt_av+=gt/(double)M;
			ht_av+=ht/(double)M;
			ct_av+=ct/(double)M;
			
			if (m%1==0)printf("%lf %d %d %lf %lf %lf %lf %lf %lf\n",q,m,cnum,el2,gl2,hl2,lcl2,gcratio,lcratio);
			fprintf(fq,"%lf %d %d %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n",q,m,cnum,el2,gl2,hl2,lcl2,gcratio,lcratio,et,gt,ht,ct);
			
			if (nb){
				for (i=0;i<N;i++){
					if (nb[i])free(nb[i]);
				}
				free(nb);
			}
			
			if (rule){
				for (i=0;i<N;i++){
					if (rule[i])free(rule[i]);
				}
				free(rule);
			}
		}//end of m-loop
		
		fprintf(fp,"%lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n",q,el2_av,gl2_av,hl2_av,lcl2_av,gcratio_av,lcratio_av,et_av,gt_av,ht_av,ct_av);
	}//end of q-loop
	
	fclose(fp);
	fclose(fq);
	
	return 0;
}

