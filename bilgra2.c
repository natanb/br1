/***********************************************************
/                                                          /
/  Programma per l'interpolazione di una griglia regolare  /
/  di punti mediante SPLINES BILINEARI                     /
/                                                          /
***********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define MMAX 2720000
#define NMAX 1905160
#define NBMAX 8
#define NNBMAX 7*NMAX
#define max(a , b) ( (a) > (b) ? (a) : (b) )
#define min(a , b) ( (a) < (b) ? (a) : (b) )

/* PROTOTYPE */

void  algite(double     *c , double   *t , int   *jc , int   *lc ,
             double     *d , int     *ic , int   *kc , int   *iw , int *ir,
             double rermax , int  itemax , int   ite , int  iter ,
             int     ifine , double   *r , double *p , double *w , int n);

void pmati(double *c, int *jc, int *lc, double *p, double *w, int n);
void classi(int *inform, int *kont, int *ipunt, int m, int n);
void inveri(double *t, int *jc, int *lc, int *ic, int *kc,
            double *w, int *iw, double *r, int *ir, double *q, int n);
void solsi(double *t, int *jc, int *lc, double *d, int n);
void cnjgrd(double *c, double *t, int *jc,int *lc,
            double *d, double rermax, int itemax,int ite,
            double *r, double *p, double *w, int n);
void tcholi(double *t, int *jc, int *lc, double *w, int n);

int offsj3(int ii,int jj,int nix);
void norsym(int *lc, int *jc, int nix, int niy, int nit);

double fmati1 (double *c, int *jc, int *lc, double *w, int n);
double fmati2 (double *c, int *ic, int *jc, int *kc, double *w, int n);

static double t[NNBMAX],c[NNBMAX];
static int jc[NNBMAX],lc[NMAX],ivive[NMAX],nng[NMAX];
static double d[NMAX],r[NMAX],u[NMAX],w[NMAX];
static int    ic[NNBMAX],kc[NMAX],iw[NMAX],ir[NMAX];

struct equaz {
          int    ia[4]; /* indici delle incognite matrice disegno */
          double a[4]; /* coeff. delle incognite matrice disegno */
          float  ob[4]; /* osservazioni x,y,z,p */
          };

/* FUNZIONE MAIN   */

main()
{
FILE *ftmpobs;
FILE *ftmpequ;
char   nf[65];
double a[4];
int    ia[4];
int    ni[5];
float  obs[4];
struct equaz  eq;

/*variabili di servizio*/
int    i,j,k,m,n,nx,ny,ncx,ncy,nit,ntc,nix,niy,ncmn,nc,ii,is,js,jss,jj,j2,jcs;
int    nii,nif,kk,i1,i2,n1,n2,n3,n4,ite,nn;
int    idd,iss,iaa,ibb;
double stp,x0,y0,x1,y1,xx,yy,zz,pp,srvx,srvy,ez,sz,s,aa,pd1,pd2,sqm;
FILE   *fdi;
float *fd;
int   *ina;
double *coef;

scanf("%lf %lf %lf %lf %lf %lf %lf",&stp, &x0, &y0, &x1, &y1, &pd1, &pd2);
scanf("%s",nf);

/*apertura file dati e file temporanei */

if( (fdi=fopen(nf,"r")) == NULL )
 {
   printf("File open %s failed \n", nf);
   exit(1);
  }

if( (ftmpobs=fopen("tmpobs","w+b")) == NULL )
{
   printf("TmpFile for observations could be not created\n");
   exit(1);
}
if( (ftmpequ=fopen("tmpequ","w+b")) == NULL )
{
   printf("TmpFile for equation (design matrix, tn...) could be not created\n");
   exit(1);
}


for( i = 1; i < NMAX; i++) {
    nng[i]=0;
}

ez=0.;
sz=0.;
i=0;
kk=0;

    while(1)
    {
    kk++;
         fscanf(fdi,"%f %f %f %f",&obs[0],&obs[1],&obs[2],&obs[3]);
         if(feof(fdi)) break;
         xx=obs[0];
         yy=obs[1];
         zz=obs[3];
	 obs[2]=zz;
         obs[3]=1.;
//	 obs[3]=pd2/(obs[3]*obs[3]);
         if( (xx > x0) && (xx < x1) && (yy > y0) && (yy < y1) ) {
            i++;
            ez+=zz;
            sz+=zz*zz;
         fwrite(obs,sizeof(float),4,ftmpobs);
         }
    }

   m=i;
   fclose (fdi);
   rewind(ftmpobs);

          ez=ez/m;
	  printf("No. Data = %d\n",kk);
          printf("M(z) = %lf \n",ez);
          printf("SD(z) = %lf \n",sqrt((sz/m-ez*ez)*m/(m-1)) );
//          printf("SD(z) = %lf \n",(sz/m-ez*ez)*m/(m-1) );

/*          if(m > MMAX) {
            printf("numero osservazioni eccedente il dimensionamento\n");
            exit(1);
          }
*/
	  x0=x0-stp/10.;
	  y0=y0-stp/10.;
	  x1=x1+stp/10.;
	  y1=y1+stp/10.;
          ncx=(int)((x1-x0)/stp);
          ncy=(int)((y1-y0)/stp);

          ntc=ncx*ncy;
          nix=ncx+1;
          niy=ncy+1;
          nit=nix*niy;
          printf("ncx = %d  ncy = %d  ntc = %d  nit= %d\n", ncx,ncy,ntc,(ncx+1)*(ncy+1));
          printf("step = %lf xmin = %lf  ymin = %lf\n",stp,x0,y0);
          printf("           xmax = %lf  ymax = %lf\n",x1,y1);
          printf("n. oss = %d \n",m);

          for( i = 1; i<=nit; i++)
          {
            lc[i]=0;
            d[i]=0.;
            ivive[i]=0;
          }
          lc[nit+1]=0;

/*
/   Costruzione della matrice normale simbolica in forma compatta
/   si considerano anche gli effetti di bordo.
*/
          norsym(lc,jc,nix,niy,nit);

          ncmn=lc[nit+1]-1;
          printf("Dimensione della normale compatta = %d\n",ncmn);

          for( i = 1; i<= ncmn; i++)
            c[i]=0.;

          for( i = 1; i<= m; i++)
          {
          fread(obs,sizeof(float),4,ftmpobs);
            xx=obs[0];
            yy=obs[1];
            zz=obs[2];
            pp=obs[3];
            xx=xx-x0;
            yy=yy-y0;
            nx=(int)(xx/stp)+1;
            ny=(int)(yy/stp)+1;
            nc=(ny-1)*ncx+nx;
            nng[nc]=nng[nc]+1;
            srvx=(xx-(nx-1)*stp)/stp;
            srvy=(yy-(ny-1)*stp)/stp;
            ia[0]=(ny-1)*nix+nx;
            ia[1]=ia[0]+1;
            ia[2]=ia[0]+nix;
            ia[3]=ia[2]+1;
            a[0]=1.-srvx-srvy+srvx*srvy;
            a[1]=srvx-srvx*srvy;
            a[2]=srvy-srvx*srvy;
            a[3]=srvx*srvy;
/*Su file temporaneo a[0...3] e ia[0...3] */
            for(k=0; k<=3; k++)
            {
                eq.ia[k]=ia[k];
                eq.a[k]= a[k];
                eq.ob[k]=obs[k];
            }
            fwrite(&eq,sizeof(struct equaz),1,ftmpequ);

            ivive[ia[0]]+=1;
            ivive[ia[1]]+=1;
            ivive[ia[2]]+=1;
            ivive[ia[3]]+=1;
            for( is = 0; is <= 3; is++)
            {
                ii=ia[is];
                aa=a[is];
                for( js = is; js <= 3; js++)
                {
                    jj=ia[js];
                    j2=lc[ia[js]+1]-1;
                    /* jcs=distanza dell'elemento c(ii,jj) da c(jj,jj) nella
                    /  normale compatta
                    */
                    jcs=offsj3(ii,jj,nix);

                    jcs=j2-jcs;
                    if(jcs == 0)
                    {
                       printf("errore %d %d oss = %d %d %d\n",j2,jcs,i,ii,jj);
                       printf("errore %lf %lf %lf %lf\n",xx,yy,zz,pp);
                       exit(1);
                    }
                    else
                    {
                     c[jcs]=c[jcs]+aa*a[js]*pp;
                    }

                 }
                 d[ii]=d[ii]-aa*zz*pp;
            }
          }

          fclose(ftmpobs);

	  idd=0;
          for(i = 1; i<=nit; i++)
           if(ivive[i] !=0) idd++;
 
          printf("nit_eff = %d \n",idd); 

/* Vincoli sul gradiente   jj >= ii*/

          for(i = 1; i<=nit; i++)
          {
            nx=(i-1)%nix + 1;
            ny=(i-nx)/nix + 1;
            idd=min(nix-nx,1);
            iss=min(nx-1,1);
            ibb=min(ny-1,1);
            iaa=min(niy-ny,1);
/* vincolo del gradiente nella direzione x */
            ii=i-iss;
            jj=i+idd;
            j2=lc[ii+1]-1;
            c[j2]=c[j2]+pd1;
            j2=lc[jj+1]-1;
            c[j2]=c[j2]+pd1;
            jcs=offsj3(ii,jj,nix);
            j2=lc[jj+1]-1-jcs;
            c[j2]=c[j2]-pd1;


/* vincolo del gradiente nella direzione y  */
            ii=i-ibb*nix;
            jj=i+iaa*nix;
            j2=lc[ii+1]-1;
            c[j2]=c[j2]+pd1;
            j2=lc[jj+1]-1;
            c[j2]=c[j2]+pd1;
            jcs=offsj3(ii,jj,nix);
            j2=lc[jj+1]-1-jcs;
            c[j2]=c[j2]-pd1;
           }

          n=nit;
          s=fmati1(c,jc,lc,w,n);
/* Soluzione con ICCG */
          algite(c,t,jc,lc,d,ic,kc,iw,ir,1.e-8,100,ite,1,1,r,u,w,n);
          aa=fmati2(c,ic,jc,kc,w,n);
	  printf("n.cond.=  %lf %lf %lf\n",s/aa,s,aa);
/*calcola scarti residui e sigma zero */
          rewind(ftmpequ);
          sz=0.;
          for(i=1; i<=m; i++)
          {
             s=0;
             fread(&eq,sizeof(struct equaz),1,ftmpequ);

             s= eq.a[0]*d[eq.ia[0]]+   /* scarto residuo sulla i-esima eq. */
                eq.a[1]*d[eq.ia[1]]+
                eq.a[2]*d[eq.ia[2]]+
                eq.a[3]*d[eq.ia[3]];
             t[i]=s-eq.ob[2];
/*             printf("%6.1f %6.1f %4.1f %4.1lf %5.2lf\n",eq.ob[0],eq.ob[1],
                     eq.ob[2],s,t[i]);
*/
             sz=sz+t[i]*t[i]*eq.ob[3];
          }
          for(i=1; i<=n; i++)
          {
            
            nx=(i-1)%nix + 1;
            ny=(i-nx)/nix + 1;
            idd=min(nix-nx,1);
            iss=min(nx-1,1);
            ibb=min(ny-1,1);
            iaa=min(niy-ny,1);
/* scarto sull'eq. vincolo del gradiente nella direzione x */
            ii=i-iss;
            jj=i+idd;
            s=d[jj]-d[ii];
            sz=sz+s*s*pd1;
/*            t[2*i-1+m]=s;*/
/*             printf("%d %4.1f %4.1f %5.2lf\n",i,d[jj],d[ii],s);*/
 
/* scarto sell'eq. vincolo del gradiente nella direzione y  */
            ii=i-ibb*nix;
            jj=i+iaa*nix;
            s=d[jj]-d[ii];
            sz=sz+s*s*pd1;
/*            t[2*i+m]=s;*/
/*             printf("%d %4.1f %4.1f %5.2lf\n",i,d[jj],d[ii],s);*/
          }  
          sz=sz/(m+nit); /* (m-n+2n) m=eq n=inc 2n=vincoli */ 
          sqm=sqrt(sz);
          printf("sigma_zero = %lf \n",sqm);
          for(i = 1; i<=nit; i++)
          {
            nx=(i-1)%nix+1;
            ny=(i-nx)/nix+1;
            xx=x0+(nx-1)*stp;
            yy=y0+(ny-1)*stp;
            s=0.;
            if(sz*c[lc[i+1]-1] > 1.e-10) s=sqrt(sz*c[lc[i+1]-1]);
          printf("%10.6lf %10.6lf %6.1lf %6.1lf %3d\n",xx,yy,d[i],s,ivive[i]);
          }

/* calcolo sqm su scarti residui */
          rewind(ftmpequ);
	  iss=0;
          for(k=1; k<=m; k++)
          {
             aa=0.;             

             fread(&eq,sizeof(struct equaz),1,ftmpequ);
             for(is=0;is<=3;is++)
             {
		s=0.;
	        i=eq.ia[is];
                i1=kc[i]+1;
                i2=kc[i+1];
                for(jss=i1;jss<=i2;jss++)
                {
                 js=ic[jss];
                 j=jc[js];
                 w[j]=c[js];
                }
                 for(js=0;js<=3;js++)
		 {
 	          j=eq.ia[js];
                  if(j >=  i) s=s+w[j]*eq.a[js];
                 }
                 aa=aa+2.*eq.a[is]*s-w[i]*eq.a[is]*eq.a[is];
              }
             pp=1./eq.ob[3]-aa;
             aa=0.;
             sqm=sqrt(pd2/eq.ob[3]);
             if(sz*pp > 1.e-10) 
             {
               aa=sqrt(sz*pp);
               if(fabs(t[k]/aa) >= 2.57) 
               {
                 iss++; 
                 printf("%6.1f %6.1f %4.1f %4.1lf %5.2lf %5.2lf  OUT\n",eq.ob[0],
                       eq.ob[1],eq.ob[2],sqm,t[k],aa);
               }
               else 
               {
                 printf("%6.1f %6.1f %4.1f %4.1lf %5.2lf %5.2lf\n",eq.ob[0],
                  eq.ob[1],eq.ob[2],sqm,t[k],aa);
               }
              }
              else
              {
                 printf("%6.1f %6.1f %4.1f %4.1lf %5.2lf %5.2lf\n",eq.ob[0],
                  eq.ob[1],eq.ob[2],sqm,t[k],aa);
              }
            }

            printf("Outliers = %d\n",iss);
fclose(ftmpequ);

}

/**********************************************************************
*                                                                      *
* nix = numero di incognite sull'asse x (input)                        *
* niy = numero di incognite sull'asse y (input)                        *
* lc  = puntatore alle colonne della matrice normale                   *
* jc  = indicatore delle righe della matrice normale                   *
* nit = numero di incognite                                            *
* Versione con dentro coefficienti del gradiente                       *
************************************************************************/
/***********************************************************************
           i+nix-1     i+nix     i+nix-1
   i-2     i-1         i         i+1
           i-nix-1     i-nix     i-nix+1
                       i-2*nix

le inc. i-2 e i-2*nix sono coinvolte nell'eq. di pseudo-oss.
Attenzione ai BORDIIIII
************************************************************************/


        void norsym(int *lc, int *jc, int nix, int niy, int nit)
{
        int iv[9];
        int k,l,kk,niv,i,j,ib,ia,is,id;
        int nx,ny,nibs;
        niv=0;
        lc[1]=1;
        kk=1;
        for(i = 1; i <= nit; i++)
        {
          for(j = 1;  j <= 6; j++)
            iv[j]=0;

          nx=(i-1)%nix + 1;
          ny=(i-nx)/nix + 1;
          id=min(nix-nx,1);
          is=min(nx-1,1);
          ib=min(ny-1,1);
          ia=min(niy-ny,1);
          nibs=i-ib*nix-is;
          ib=ib+1;
          k=0;

          if(min(ny-1,2) == 2) {
            k=k+1;
            iv[k]=nibs+is-nix;
            }
          for(l = 1; l <= ib; l++)
            {
             if(l==ib && min(nx-1,2) == 2) {
               k=k+1;
               iv[k]=nibs-1;
               }
             for(j = 0; j <= is+id; j++)
              {
              k=k+1;
              iv[k]=nibs+j;
              }
            nibs=nibs+nix;
            }
          niv=0;

          for(l = k; l >= 1; l--)
            if(iv[l] == i) niv=l;

          if(niv == 0)
            {
            printf("errore nella costruzione simbolica");
            exit(1);
            }

          lc[i+1]=lc[i]+niv;
          for(j = 1; j <= niv; j++)
            {
            jc[kk]=iv[j];
            kk=kk+1;
            }

        }
        return;
}

/********************************************************/
int offsj3(int ii, int jj, int nix)
{
    int is,js,ib,jb,im,jm,jcs;
    is=ii-1;
    js=jj-1;
    ib=is/nix;
    jb=js/nix;
    im= is % nix + 1;
    jm= js % nix + 1;
    if(jb == ib ) return (jm-im);
    if((jb-ib) == 1) return ( min(jm,nix-1)-im+2 + min(jm,3) - 1 );
    if((jb-ib) == 2) return ( min(2*jm+1,6)-max(jm-nix+1,0) );
/*    if(jb-ib == 2) return ( min(2*jm+1,6)-jm/nix );
L'espressione commentata e' equivalente infatti
jm/nix = max(jm-nix+1,0)!!!!!!!!!!!!!!                      */
}

/***********************************************************/

void  algite(double *c, double *t, int *jc, int *lc,
              double *d, int *ic, int *kc, int *iw, int *ir,
             double rermax, int itemax, int ite, int iter,
             int ifine, double *r, double *p, double *w, int n)
{
int     i,j,nn;


   if( iter == 1 )
   {
      nn=lc[n+1]-1;
      for( i = 1; i <= nn; i++)
         t[i]=c[i];

      tcholi(t,jc,lc,w,n);
   }

   cnjgrd(c,t,jc,lc,d,rermax,itemax,ite,r,p,w,n);

   if(ifine != 1) return;
   inveri(t,jc,lc,ic,kc,w,iw,r,ir,p,n);
   nn=lc[n+1]-1;
   for(i=1; i<=nn;i++)
   c[i]=t[i];
   return;
}

/*************************          tcholi          *************************/
void tcholi(double *t, int *jc, int *lc, double *w, int n)
{
int     i,i1,i2;
int     j,j1,j2,js;
int     k,ks;
double  s;

 for(i = 1; i <= n; i++)
  w[i]=0.;

for(i = 1; i <= n; i++)
 {
   i1=lc[i];
   i2=lc[i+1]-1;
   for(js = i1; js <= i2; js++)
   {
     j=jc[js];
     j1=lc[j];
     j2=lc[j+1]-2;
     s=t[js];
     if(j1 <= j2)
     {
      for(ks = j1; ks <= j2; ks++)
      {
        k=jc[ks];
        s=s-t[ks]*w[k];
      }
     }
     t[js]=s;
     j2=j2+1;
     w[j]=s/t[j2];
   }
   for(js = i1; js <= i2; js++)
   {
     j=jc[js];
     w[j]=0.;
   }
 }
 return;
}

/*************************          cnjgrd          *************************/
void cnjgrd(double *c, double     *t, int    *jc, int *lc,
            double *d, double rermax, int itemax, int ite,
            double *r, double     *p, double  *w, int n)
{
int     i;
double  s,rtr,rmax,dmax;
double  ptcp,alfa,beta;

          for(i = 1; i <= n; i++)
          {
             r[i]=-d[i];
             p[i]=r[i];
          }
          solsi(t,jc,lc,p,n);
          s=0.0;
          for(i = 1; i <= n; i++)
          {
             s=s+r[i]*p[i];
             d[i]=0.0;
          }
          rtr=s;
          ite=0;
loop:     rmax=0.0;
          dmax=0.0;
          for(i = 1; i <= n; i++)
          {
            rmax=max(rmax,fabs(r[i]));
            dmax=max(dmax,fabs(d[i]));
          }
          printf("ite = %d\n",ite);
          if(rmax <= (dmax*rermax) || ite >= itemax) return;
          ite=ite+1;
          pmati(c,jc,lc,p,w,n);
          s=0.0;
          for(i = 1; i <= n; i++)
             s=s+p[i]*w[i];

          ptcp=s;
          alfa=rtr/ptcp;
          for(i = 1; i <= n; i++)
          {
             d[i]=d[i]+alfa*p[i];
             r[i]=r[i]-alfa*w[i];
             w[i]=r[i];
          }
          solsi(t,jc,lc,w,n);

          s=0.0;
          for(i = 1; i <= n; i++)
             s=s+r[i]*w[i];

          beta=s/rtr;
          rtr=s;
          for(i = 1; i <= n; i++)
            p[i]=w[i]+beta*p[i];

          goto loop;
}

/*************************          solsi            *************************/
void solsi(double *t, int *jc, int *lc, double *d, int n)
{
int i,i1,i2,is;
int j,js;
double s,ds;

 for(i = 1; i <= n; i++)
 {
   i1=lc[i];
   i2=lc[i+1]-2;
   s=d[i];
   if(i1 <= i2)
   {
    for(js = i1; js <= i2; js++)
    {
      j=jc[js];
      s=s-t[js]*d[j];
    }
   }
   i2=i2+1;
   d[i]=s/t[i2];
 }
 for(i = 1; i <= n; i++)
 {
   is=lc[i+1]-1;
   d[i]=d[i]*t[is];
 }
 for(i = n; i >= 1; i--)
 {
   i1=lc[i];
   i2=lc[i+1]-1;
   d[i]=d[i]/t[i2];
   i2=i2-1;
   if(i1 <= i2)
   {
    ds=d[i];
    for(js = i2; js >= i1; js--)
    {
     j=jc[js];
     d[j]=d[j]-t[js]*ds;
    }
   }
  }
 return;
}

/*************************          pmati            *************************/
void pmati(double *c, int *jc, int *lc, double *p, double *w, int n)
{
int     i,i1,i2;
int     j,js;
double  s,ps;

 for(i = 1; i <= n; i++)
 {
   i1=lc[i];
   i2=lc[i+1]-1;
   s=0.0;
   for(js = i1; js <= i2; js++)
   {
     j=jc[js];
     s=s+c[js]*p[j];
   }
   w[i]=s;
 }
 for(i = n; i >= 1; i--)
 {
   i1=lc[i];
   i2=lc[i+1]-2;
   if(i1 <= i2)
   {
    ps=p[i];
    for(js = i2; js >= i1; js--)
    {
      j=jc[js];
      w[j]=w[j]+c[js]*ps;
    }
   }
 }
 return;
}

/*************************          inveri          *************************/
void inveri(double *t, int *jc, int *lc, int *ic, int *kc,
            double *w, int *iw, double *r, int *ir, double *q, int n)
{
int     i,i1,i2,is;
int     j,j1,j2,jj,js;
int     nn,l,ll,ls,lss,k,ks,kk;
double  qq,tt,rr,s;


 nn=lc[n+1]-1;
 t[nn]=1.0/t[nn];
 classi(jc,kc,ic,nn,n);

 for(i = 1; i <= n; i++)
 {
   w[i]=0.0;
   r[i]=0.0;
   q[i]=0.0;
   i1=lc[i];
   i2=lc[i+1]-1;
   for(js = i1; js <= i2; js++)
     jc[js]=i;
 }

 for(i = n-1; i >= 1; i--)
 {
   i1=kc[i]+1;
   i2=kc[i+1];
   is=ic[i1];
   tt=t[is];
   qq=0.0;
   i1=i1+1;
   if( i1 <= i2)
   {
    kk=0;
    for(ks = i1; ks <= i2; ks++)
    {
      js=ic[ks];
      j=jc[js];
      kk=kk+1;
      ir[kk]=j;
      r[j]=t[js];
    }
    k=0;
    jj=0;
    for(ks = i1; ks <= i2; ks++)
    {
      js=ic[ks];
      k=k+1;
      j=ir[k];
      j1=kc[j]+1;
      j2=kc[j+1];
      ll=0;
      for(lss = j1; lss <= j2; lss++)
      {
        js=ic[ks];
        ls=ic[lss];
        l=jc[ls];
        ll=ll+1;
        iw[ll]=l;
        w[l]=t[ls];
      }
      jj=max(jj,l);
      rr=r[j];
      s=0.0;
      for(ls = 1; ls <= ll; ls++)
      {
        js=ic[ks];
        l=iw[ls];
        s=s+r[l]*w[l];
        q[l]=q[l]+rr*w[l];
      }
      js=ic[ks];
      t[js]=-(q[j]-rr*w[j]+s)/tt;
      qq=qq+rr*t[js];
      for(ls = 1; ls <= ll; ls++)
      {
        l=iw[ls];
        w[l]=0.0;
      }
    }

    for(ks = 1; ks <= kk; ks++)
    {
      j=ir[ks];
      r[j]=0.0;
    }
    for(j = i; j <= jj; j++)
      q[j]=0.0;

   }
   t[is]=(1.0-qq)/tt;
 }
 return;
}

/******************************     classi **************************/
void classi(int *inform, int *kont, int *ipunt, int m, int n)
{
 int     i,k,k1,k2;
 for(i = 1; i <= n; i++)
   kont[i+1]=0;

 kont[1]=0;
 for(i = 1; i <= m; i++)
 {
   k=inform[i];
   if(k != 0)
    kont[k+1]=kont[k+1]+1;
 }

 for(i = 1; i <= n; i++)
   kont[i+1]=kont[i]+kont[i+1];

 for(i = 1; i <= m; i++)
 {
   k1=inform[i];
   if(k1 != 0)
   {
    kont[k1]=kont[k1]+1;
    k2=kont[k1];
    ipunt[k2]=i;
   }
 }

 for(i = n; i >= 1; i--)
   kont[i+1]=kont[i];

 kont[1]=0;
 return;
}
/*************************          fmati          *************************/
double fmati1 (double *c, int *jc, int *lc, double *w, int n)
{
  int i, is, j, js, ini, ifi;
  double wmax;

  for (i = 1; i <= n; i++)
    w[i] = 0.;

  w[1] = c[1];

  for (i = 2; i <= n; i++)
    {
      is = lc[i + 1] - 1;
      w[i] = w[i] + c[is];
      ini = lc[i];
      ifi = is - 1;
      if (ifi >= ini)
	for (js = ini; js <= ifi; js++)
	  {
	    j = jc[js];
	    w[i] = w[i] + fabs(c[js]);
	    w[j] = w[j] + fabs(c[js]);
	  }
    }


  wmax = 0.;
  for (i = 1; i <= n; i++)
    if (wmax <= w[i])
      wmax = w[i];
  return (wmax);
}

double fmati2 (double *c, int *ic, int *jc, int *kc, double *w, int n)
{
  int i, is, iss, j, js, jss, ini, ifi, ns, nns;
  double wmax;

  for (i = 1; i <= n; i++)
    w[i] = 0.;

  for (i = 1; i < n; i++)
    {
      iss = kc[i] + 1;
      is = ic[iss];
      w[i] = w[i] + c[is];
      ini = iss + 1;
      ifi = kc[i + 1];
      if (ifi >= ini)
	for (jss = ini; jss <= ifi; jss++)
	  {
	    js = ic[jss];
	    j = jc[js];
	    w[i] = w[i] + fabs(c[js]);
	    w[j] = w[j] + fabs(c[js]);
	  }
    }
  nns = kc[n] + 1;
  ns = ic[nns];
  w[n] = w[n] + c[ns];
  wmax = 0.;

    for (i = 1; i <= n; i++)
    if (wmax < w[i])
      wmax = w[i];
  return (wmax);
}
