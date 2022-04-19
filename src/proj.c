#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sched.h>
#include <omp.h>
#include "allocate.h"
#include "data.h"
#include "prepro.h"
#include "proj.h"
#include "mpi.h"
#include <string.h>

void createSinogram(struct Sinogram *sinogram)
{
	sinogram->sino = (ENTRY *)get_spc((sinogram->geom_info.Nr)*(sinogram->geom_info.Nc)*(sinogram->geom_info.Nv), sizeof(ENTRY));
}

void freeSinogram(struct Sinogram *sinogram)
{

	fprintf(stdout,"inside freeSinogram \n");
	fflush(stdout);

	if(sinogram->sino != NULL)
		free(sinogram->sino);

	fprintf(stdout,"after free sino \n");
	fflush(stdout);

	if(sinogram->counts != NULL)
        	free(sinogram->counts);
	if(sinogram->dose != NULL)
		free(sinogram->dose);         


	fprintf(stdout,"after free counts \n");
	fflush(stdout);
 
        if(sinogram->offset != NULL)
		free(sinogram->offset);
        if(sinogram->D != NULL)
        	free(sinogram->D); 
}

void createImage(struct Image *image)
{
	image->img = (ENTRY *)get_spc((image->img_info.Nx)*(image->img_info.Ny)*(image->img_info.Nz), sizeof(ENTRY));
}

void freeImage(struct Image *image)
{
	free(image->img);
}

void fillGeomInfo(struct GeomInfo *geom_info)	/* fill in the intermediate variables */
{
	geom_info->alphac0 = -(geom_info->Nc-1.0)*(geom_info->Del_alphac)/2.0 + geom_info->del_alphac;
	geom_info->Del_zs = (geom_info->u)*(geom_info->Del_dr)*(fabs(geom_info->Del_beta))/(2.0*PI);
	geom_info->detc = (geom_info->Nc)*(geom_info->Del_alphac);
	geom_info->detr = (geom_info->Nr)*(geom_info->Del_dr);
	geom_info->half_detr = (geom_info->Nr-1.0)*(geom_info->Del_dr)/2.0 + geom_info->offset_dr;   /*it's actually half detector minus half a channel */
	geom_info->cone_zbuffer= geom_info->detr*(geom_info->r_si+geom_info->fov/2.0)/(2.0*geom_info->r_sd);


	//if(myid==0){
	//	fprintf(stdout,"alphac0 %f Del_zs %f detc %f detr %f half_detr %f \n",geom_info->alphac0,geom_info->Del_zs,geom_info->detc,geom_info->detr,geom_info->half_detr);
	//}

}

void fillImgInfo(struct ImgInfo *img_info)	/* fill in the intermediate variables */
{
	img_info->x0 = img_info->xc - (img_info->Del_xy)*(img_info->Nx-1)/2.0;
	img_info->y0 = img_info->yc - (img_info->Del_xy)*(img_info->Ny-1)/2.0;
	img_info->z0 = img_info->zc - (img_info->Del_z)*(img_info->Nz-1)/2.0;

	//fprintf(stdout,"x0 %f y0 %f z0 %f \n",img_info->x0,img_info->y0,img_info->z0);
}

/* sjk: check information for consistency */
void checkInfo(struct GeomInfo *geom_info, struct ImgInfo *img_info)
{
        float x,y,dist;

        /* find pixel distance furthest from iso (just check corners) */
        x = img_info->x0 - img_info->Del_xy/2.0;
        y = img_info->y0 - img_info->Del_xy/2.0;
        dist = sqrt(x*x + y*y);

        x = img_info->x0 - img_info->Del_xy/2.0 + img_info->Nx;
        y = img_info->y0 - img_info->Del_xy/2.0;
        dist = max(dist,sqrt(x*x+y*y));

        x = img_info->x0 - img_info->Del_xy/2.0;
        y = img_info->y0 - img_info->Del_xy/2.0 + img_info->Ny;
        dist = max(dist,sqrt(x*x+y*y));

        x = img_info->x0 - img_info->Del_xy/2.0 + img_info->Nx;
        y = img_info->y0 - img_info->Del_xy/2.0 + img_info->Ny;
        dist = max(dist,sqrt(x*x+y*y));

        if(dist>geom_info->r_si) {
                fprintf(stdout,"ERROR: You have pixels outside the source radius.");
                exit(-1);
        }
}


void createSourceLocInfo(
	struct SourceLocInfo *source_loc_info,
	struct GeomInfo *geom_info)
{
	source_loc_info->beta = (ENTRY *)get_spc(geom_info->Nv, sizeof(ENTRY));
	source_loc_info->xs = (ENTRY *)get_spc(geom_info->Nv, sizeof(ENTRY));
	source_loc_info->ys = (ENTRY *)get_spc(geom_info->Nv, sizeof(ENTRY));
	source_loc_info->zs = (ENTRY *)get_spc(geom_info->Nv, sizeof(ENTRY));
}

void compSourceLocInfo(
	struct SourceLocInfo *source_loc_info,
	struct GeomInfo *geom_info,int myid,int total_nodes)
{
	int iv;		/* view index */
	int i,j;

	FILE* fp;
   
	fprintf(stdout,"reach compSourceLocInfo\n");
	fflush(stdout);
 
    	if ((fp = fopen(geom_info->viewAnglesList, "r")) == NULL)
    	{
        fprintf(stderr, "ERROR! Can't open file containing list of view angles %s.\n", geom_info->viewAnglesList);
        exit(-1);
    	}

    	for(j=0;j< (myid/(geom_info->num_focal_spots*geom_info->num_sources));j++){
		fscanf(fp, "%*[^\n]\n");
    	}

    
    	for(i=0;i<geom_info->Nv;i++)
    	{
        	if(fscanf(fp,"%f\n",&(source_loc_info->beta[i])) == 0)
       		{
         		fprintf(stderr, "ERROR! List of view angles in file %s terminated early.\n", geom_info->viewAnglesList);
         		exit(-1);
       		}
		for(j=0;j< ((total_nodes/(geom_info->num_focal_spots*geom_info->num_sources))-1);j++){
			fscanf(fp, "%*[^\n]\n");
		}
    	}



	if ((fp = fopen(geom_info->zPositionList, "r")) == NULL)
    	{
        	fprintf(stderr, "ERROR! Can't open file containing list of source z position %s.\n", geom_info->zPositionList);
        	exit(-1);
    	}

    	for(j=0;j< (myid/(geom_info->num_focal_spots*geom_info->num_sources));j++){
		fscanf(fp, "%*[^\n]\n");
    	}


    
    	for(i=0;i<geom_info->Nv;i++)
    	{
        	if(fscanf(fp,"%f\n",&(source_loc_info->zs[i])) == 0)
       		{
         		fprintf(stderr, "ERROR! List of source z position in file %s terminated early.\n", geom_info->zPositionList);
         		exit(-1);
       		}
		for(j=0;j< ((total_nodes/(geom_info->num_focal_spots*geom_info->num_sources))-1);j++){
			fscanf(fp, "%*[^\n]\n");
		}
    	}
    
    	fclose(fp);



	//if(myid==3){
	//for (iv =0; iv< geom_info->Nv; iv++){
	//		fprintf(stdout,"iv %d zs %f \n",iv,source_loc_info->zs[iv]);
	//}
	//}



	for (iv = 0; iv < geom_info->Nv; iv++)
	{
		source_loc_info->xs[iv] = (geom_info->r_si)*cos(source_loc_info->beta[iv]);
		source_loc_info->ys[iv] = (geom_info->r_si)*sin(source_loc_info->beta[iv]);
		
	}


	fprintf(stdout,"finish compute source loc \n");
	fflush(stdout);


/*	
	for (iv=0;iv<geom_info->Nv;iv++){
		if(myid==2){
			fprintf(stdout,"iv %d zs %f beta %f \n",iv,source_loc_info->zs[iv],source_loc_info->beta[iv]);
		}
	}	
	fflush(stdout);
*/



}

void freeSourceLocInfo(struct SourceLocInfo *source_loc_info)
{
	free(source_loc_info->beta);
	free(source_loc_info->xs);
	free(source_loc_info->ys);
	free(source_loc_info->zs);
}




void createViewXYInfo_stored(
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info, struct ImgInfo * img_info,char **recon_mask)
{

	view_xy_info->ic_start = (CHANNEL ***)malloc(sizeof(CHANNEL**) *img_info->Nx);
	view_xy_info->ic_num = (PROCHANNEL ***)malloc( sizeof(PROCHANNEL**) *img_info->Nx);
	view_xy_info->Mag = (unsigned char ***)malloc( sizeof(unsigned char**)*img_info->Nx);
	//view_xy_info->Wr = (ENTRY ***)malloc(sizeof(ENTRY**)*img_info->Nx);
	view_xy_info->B = (unsigned char ****)malloc(sizeof(unsigned char***)*img_info->Nx);

	view_xy_info->Mag_MaxPointer = (ENTRY *)malloc(sizeof(ENTRY)*img_info->Nx *img_info->Ny);
	view_xy_info->B_MaxPointer = (ENTRY *)malloc(sizeof(ENTRY)*img_info->Nx *img_info->Ny);
	//fprintf(stdout,"after allocating 1 \n");
	//fflush(stdout);

	for(int jx=0;jx<img_info->Nx;jx++){

		view_xy_info->ic_start[jx] = (CHANNEL **)malloc(sizeof(CHANNEL*) *img_info->Ny);
		view_xy_info->ic_num[jx] = (PROCHANNEL **)malloc( sizeof(PROCHANNEL*) *img_info->Ny);
		view_xy_info->Mag[jx] = (unsigned char **)malloc( sizeof(unsigned char*)*img_info->Ny);
		//view_xy_info->Wr[jx] = (ENTRY **)malloc(sizeof(ENTRY*)*img_info->Ny);
		view_xy_info->B[jx] = (unsigned char ***)malloc(sizeof(unsigned char**)*img_info->Ny);

		//fprintf(stdout,"after allocating jx %d \n",jx);
		//fflush(stdout);

		for (int jy=0;jy<img_info->Ny;jy++){
		if(recon_mask[jx][jy]){
			view_xy_info->ic_start[jx][jy]=(CHANNEL *)malloc(sizeof(CHANNEL)*geom_info->Nv);
			view_xy_info->ic_num[jx][jy]=(PROCHANNEL *)malloc(sizeof(PROCHANNEL)*geom_info->Nv);
			view_xy_info->Mag[jx][jy] = (unsigned char *)malloc(sizeof(unsigned char)*geom_info->Nv);
			//view_xy_info->Wr[jx][jy] = (ENTRY *)malloc(sizeof(ENTRY)*geom_info->Nv);
			view_xy_info->B[jx][jy] = (unsigned char **)malloc(sizeof(unsigned char*)*geom_info->Nv);
		}
		}
	}
}


/*
void compViewXYInfo_OnTheFly(
	int jx,
	int jy,
	float x,
	float y,
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info,
	struct SourceLocInfo *source_loc_info)
{
	int iv;		
	float theta;		
	float theta_td;	
	float costh;		
	float alphaj;		
	float alphaj_td;	
	//float alpha_min;
	//float alpha_max;
	float r_sv;		
	float Wc;		
	float del_c;		
	float Bij;
	//CHANNEL ic_end;
	PROCHANNEL p;
	int ic;		






	for (iv = 0; iv < geom_info->Nv; iv++)
	{
		theta = atan2((source_loc_info->ys[iv]-y), (source_loc_info->xs[iv]-x));
		
		if (theta >= -PI/4.0)
		{
			theta_td = fmod((theta + PI/4.0), (PI/2.0)) - (PI/4.0);
		}
		else
		{
			theta_td = fmod((theta + PI/4.0), (PI/2.0)) + (PI/4.0);
		}
		costh = cos(theta_td);

		alphaj = theta - source_loc_info->beta[iv];
		alphaj_td = adjust(alphaj);

		r_sv = sqrt((source_loc_info->xs[iv]-x)*(source_loc_info->xs[iv]-x) + (source_loc_info->ys[iv]-y)*(source_loc_info->ys[iv]-y));
		
		//view_xy_info->Mag[iv] = (geom_info->r_sd)/r_sv;


		Wc = (img_info->Del_xy)*costh/r_sv;


		//alpha_min = alphaj_td - geom_info->alphac0 - (Wc - geom_info->Del_alphac)/2.0;
		//alpha_max = alphaj_td - geom_info->alphac0 + (Wc + geom_info->Del_alphac)/2.0;


		
		//if (alpha_max < 0 || alpha_min > geom_info->detc)
		//{
		//	view_xy_info->ic_num[iv] = 0;
		//}
		
		//else
		//{
			//view_xy_info->ic_start[jx*img_info->Ny+jy][iv] = (CHANNEL)max((CHANNEL)floor(alpha_min/(geom_info->Del_alphac)), 0);
			//ic_end = (CHANNEL)min((CHANNEL)floor(alpha_max/(geom_info->Del_alphac)), (CHANNEL)(geom_info->Nc-1));
			//view_xy_info->ic_num[iv] = ((PROCHANNEL)(ic_end - view_xy_info->ic_start[jx*img_info->Ny+jy][iv] + 1));
		//}

		



		view_xy_info->B[iv] = (ENTRY *)get_spc((int)(view_xy_info->ic_num[jx*img_info->Ny+jy][iv]), sizeof(ENTRY));
		
		for (p = 0; p < view_xy_info->ic_num[jx*img_info->Ny+jy][iv]; p++)
		{
			ic = (int)(view_xy_info->ic_start[jx*img_info->Ny+jy][iv] + p);

			del_c = adjust(alphaj - (ic*(geom_info->Del_alphac) + geom_info->alphac0));

			Bij = clip(0.0, ((Wc+(geom_info->Del_alphac))/2.0)-fabs(del_c), min(Wc, (geom_info->Del_alphac)));
			Bij *= ((img_info->Del_xy)/((geom_info->Del_alphac)*costh));

			view_xy_info->B[iv][p] = Bij;
		}
		
	}
}
*/




void compViewXYInfo(
	int jx,
	int jy,
	float x,
	float y,
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info,
	struct SourceLocInfo *source_loc_info)
{
	int iv;			/* view index */
	float theta;		/* slope of the ray */
	float theta_td;	/* adjusted theta, in [-PI/4, PI/4] */
	float costh;		/* cos(theta_td) */
	float alphaj;		/* voxel angle relative to ray through iso */
	float alphaj_td;	/* alphaj_td = (alphaj + PI) mod 2*PI - PI */
	float alpha_min;
	float alpha_max;
	float r_sv;		/* source to voxel */
	float Wc;		/* projection angle width */
	float del_c;		/* angle between ray through voxel center and ray through detector center */
	float Bij;
	CHANNEL ic_end;
	PROCHANNEL p;
	int ic;			/* channel index */

	ENTRY *temp_Mag = (ENTRY *)malloc(sizeof(ENTRY)*geom_info->Nv);
	
	ENTRY max_num=0;
	for (iv = 0; iv < geom_info->Nv; iv++)
	{

		r_sv = sqrt((source_loc_info->xs[iv]-x)*(source_loc_info->xs[iv]-x) + (source_loc_info->ys[iv]-y)*(source_loc_info->ys[iv]-y));
		
		temp_Mag[iv] = (geom_info->r_sd)/r_sv;

		if(temp_Mag[iv]>max_num)
			max_num=temp_Mag[iv];
	}

	view_xy_info->Mag_MaxPointer[jx*img_info->Ny+jy]=max_num;

	for (iv = 0; iv < geom_info->Nv; iv++){
       		view_xy_info->Mag[jx][jy][iv] = (unsigned char)((temp_Mag[iv])/max_num*255+0.5);
	}
	

	free((void *)temp_Mag);


	ENTRY **temp_B = (ENTRY **)malloc(sizeof(ENTRY *)*geom_info->Nv);


	for (iv = 0; iv < geom_info->Nv; iv++)
	{
		theta = atan2((source_loc_info->ys[iv]-y), (source_loc_info->xs[iv]-x));
		
		if (theta >= -PI/4.0)
		{
			theta_td = fmod((theta + PI/4.0), (PI/2.0)) - (PI/4.0);
		}
		else
		{
			theta_td = fmod((theta + PI/4.0), (PI/2.0)) + (PI/4.0);
		}
		costh = cos(theta_td);

		alphaj = theta - source_loc_info->beta[iv];
		alphaj_td = adjust(alphaj);

		//r_sv = sqrt((source_loc_info->xs[iv]-x)*(source_loc_info->xs[iv]-x) + (source_loc_info->ys[iv]-y)*(source_loc_info->ys[iv]-y));
		
		//view_xy_info->Mag[jx][jy][iv] = (geom_info->r_sd)/r_sv;

		//view_xy_info->Wr[jx][jy][iv] = (img_info->Del_z)*(view_xy_info->Mag[jx][jy][iv]);

		Wc = (img_info->Del_xy)*costh*view_xy_info->Mag[jx][jy][iv]*max_num*1.0/255/geom_info->r_sd;


		alpha_min = alphaj_td - geom_info->alphac0 - (Wc - geom_info->Del_alphac)/2.0;
		alpha_max = alphaj_td - geom_info->alphac0 + (Wc + geom_info->Del_alphac)/2.0;


		
		if (alpha_max < 0 || alpha_min > geom_info->detc)
		{
			view_xy_info->ic_num[jx][jy][iv] = 0;
		}
		
		else
		{
			view_xy_info->ic_start[jx][jy][iv] = (CHANNEL)max((CHANNEL)floor(alpha_min/(geom_info->Del_alphac)), 0);
			ic_end = (CHANNEL)min((CHANNEL)floor(alpha_max/(geom_info->Del_alphac)), (CHANNEL)(geom_info->Nc-1));
			view_xy_info->ic_num[jx][jy][iv] = ((PROCHANNEL)(ic_end - view_xy_info->ic_start[jx][jy][iv] + 1));
		}


		temp_B[iv] = (ENTRY *)malloc(sizeof(ENTRY )*((int)view_xy_info->ic_num[jx][jy][iv]));

		view_xy_info->B[jx][jy][iv] = (unsigned char *)get_spc((int)(view_xy_info->ic_num[jx][jy][iv]), sizeof(unsigned char));
	
		for (p = 0; p < view_xy_info->ic_num[jx][jy][iv]; p++)
		{
			ic = (int)(view_xy_info->ic_start[jx][jy][iv] + p);

			del_c = adjust(alphaj - (ic*(geom_info->Del_alphac) + geom_info->alphac0));

			Bij = clip(0.0, ((Wc+(geom_info->Del_alphac))/2.0)-fabs(del_c), min(Wc, (geom_info->Del_alphac)));
			Bij *= ((img_info->Del_xy)/((geom_info->Del_alphac)*costh));


			//fprintf(stdout,"jx %d jy %d iv %d Mag %f ic_start %d ic_num %d ic %d del_c %f Bij %f\n",jx,jy,iv,view_xy_info->Mag[jx][jy][iv]*max_num*1.0/255,view_xy_info->ic_start[jx][jy][iv],view_xy_info->ic_num[jx][jy][iv],ic,del_c, Bij);
			//fflush(stdout);


			temp_B[iv][p] = Bij;
		}
	}
	
	max_num=0;
	for (iv = 0; iv < geom_info->Nv; iv++)
	{
		for (p=0;p<view_xy_info->ic_num[jx][jy][iv];p++){
			if(temp_B[iv][p]>max_num)
				max_num=temp_B[iv][p];
		}	
	}
	
	view_xy_info->B_MaxPointer[jx*img_info->Ny+jy]=max_num;

	for (iv = 0; iv < geom_info->Nv; iv++)
	{
		for (p=0;p<view_xy_info->ic_num[jx][jy][iv];p++){
			view_xy_info->B[jx][jy][iv][p] =  (unsigned char)((temp_B[iv][p])/max_num*255+0.5);

		}	
	}


	for (iv = 0; iv < geom_info->Nv; iv++)
	{
		free((void *)temp_B[iv]);
	}
	free((void **)temp_B);
}




void freeViewXYInfo_stored(
	struct ViewXYInfo *view_xy_info,struct GeomInfo *geom_info,struct ImgInfo *img_info,char ** recon_mask)
{
	for(int jx=0;jx<img_info->Nx;jx++){
		for(int jy=0;jy< img_info->Ny ; jy++){
			if(recon_mask[jx][jy]){
				free(view_xy_info->ic_start[jx][jy]);
				free(view_xy_info->ic_num[jx][jy]);
				free(view_xy_info->Mag[jx][jy]);
				//free(view_xy_info->Wr[jx][jy]);
				for (int iv = 0; iv < geom_info->Nv; iv++)
				{
					free(view_xy_info->B[jx][jy][iv]);
				}
				free(view_xy_info->B[jx][jy]);
			}
		}
		free(view_xy_info->ic_start[jx]);
		free(view_xy_info->ic_num[jx]);
		free(view_xy_info->Mag[jx]);
		free(view_xy_info->B[jx]);
	}
	free(view_xy_info->ic_start);
	free(view_xy_info->ic_num);
	free(view_xy_info->Mag);
	free(view_xy_info->B);
	free(view_xy_info->Mag_MaxPointer);
	free(view_xy_info->B_MaxPointer);
}




void compViewXYZInfo(
	int jx,
	int jy,
	int jz,
	float z,
	struct ViewXYZInfo *view_xyz_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info,
	struct SourceLocInfo *source_loc_info,
	struct ViewXYInfo *view_xy_info)
{
	int iv;		/* view index */
	float d_1,d_2;
	float d_min;
	float d_max;
	ROW ir_end;
	int iv_end;  /* sjk */
	float inverse_num=1.0/255;
	int myid;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);



  	 /* sjk: this block is replaced with the lines following */
  	//{
		d_1 = (geom_info->Nr*geom_info->Del_dr/2 -geom_info->offset_dr)*(geom_info->r_si+geom_info->fov/2.0)/geom_info->r_sd;
		d_2 = (geom_info->Nr*geom_info->Del_dr/2 +geom_info->offset_dr)*(geom_info->r_si+geom_info->fov/2.0)/geom_info->r_sd;
	
		if ((source_loc_info->zs[0] <= source_loc_info->zs[geom_info->Nv-1])&&(z < (source_loc_info->zs[0]-d_1) || z > (source_loc_info->zs[geom_info->Nv-1]+d_2)))
		{
			view_xyz_info->iv_num = 0;
		}
		else if ((source_loc_info->zs[0] > source_loc_info->zs[geom_info->Nv-1])&&(z < (source_loc_info->zs[geom_info->Nv-1]-d_1) || z > (source_loc_info->zs[0]+d_2)))
		{
			view_xyz_info->iv_num = 0;
		}

		else
		{
			for (iv = 0; iv < geom_info->Nv; iv++)
			{
				if (z >= (source_loc_info->zs[iv]-d_1) && z <= (source_loc_info->zs[iv]+d_2))
				{
					view_xyz_info->iv_start = iv;
					break;
				}
			}
			view_xyz_info->iv_num = 0;

			for (iv = 0; iv < geom_info->Nv; iv++)
			{
				if (z >= (source_loc_info->zs[iv]-d_1) && z <= (source_loc_info->zs[iv]+d_2))
				{
					view_xyz_info->iv_num++;
				}
			}
		}
  	//}



	/* sjk: this block replaces the above, finding the iv range in closed form */
	/*
        view_xyz_info->iv_start=max(0,(int)floor((z-img_info->Del_z/2.0-geom_info->cone_zbuffer - geom_info->zs_0)/geom_info->Del_zs));
        iv_end =   min(geom_info->Nv-1,(int)ceil((z+img_info->Del_z/2.0+geom_info->cone_zbuffer - geom_info->zs_0)/geom_info->Del_zs));
        if((iv_end<0) || (view_xyz_info->iv_start > (geom_info->Nv-1)))
                view_xyz_info->iv_num=0;
        else
                view_xyz_info->iv_num= iv_end - view_xyz_info->iv_start + 1;
	*/

	/*
	if(temp_iv_num != view_xyz_info->iv_num || temp_iv_start != view_xyz_info->iv_start){
		fprintf(stdout,"z is %f temp_iv_num %i kisner_iv_num %i temp_iv_start %d kisner_iv_start %d \n",z,temp_iv_num,view_xyz_info->iv_num,temp_iv_start,view_xyz_info->iv_start);
		fflush(stdout);
	}
	*/

	/* sjk: moved this block down so that we don't have to go through all views */
	for (iv = view_xyz_info->iv_start; iv <= (view_xyz_info->iv_start + view_xyz_info->iv_num-1); iv++)   /* sjk */
	{
		d_min = (geom_info->half_detr+geom_info->Del_dr/2.0-view_xy_info->Mag[jx][jy][iv]*view_xy_info->Mag_MaxPointer[jx*img_info->Ny+jy]*inverse_num*(z-source_loc_info->zs[iv]+img_info->Del_z/2.0));
		d_max = (geom_info->half_detr+geom_info->Del_dr/2.0-view_xy_info->Mag[jx][jy][iv]*view_xy_info->Mag_MaxPointer[jx*img_info->Ny+jy]*inverse_num*(z-source_loc_info->zs[iv]-img_info->Del_z/2.0));

		if (d_max < 0 || d_min > geom_info->detr)
		{
			view_xyz_info->ir_num[iv] = 0;
		}
		else
		{
			view_xyz_info->ir_start[iv] = (ROW)max((ROW)floor(d_min/(geom_info->Del_dr)), 0);
			ir_end = (ROW)min((ROW)floor(d_max/(geom_info->Del_dr)), (ROW)(geom_info->Nr-1));
			view_xyz_info->ir_num[iv] = ((PROROW)(ir_end - view_xyz_info->ir_start[iv] + 1));
		}
	}
}



void createViewXYZInfo(
	struct ViewXYZInfo *view_xyz_info,
	struct GeomInfo *geom_info)
{
	view_xyz_info->ir_start = (ROW *)get_spc(geom_info->Nv, sizeof(ROW));
	view_xyz_info->ir_num = (PROROW *)get_spc(geom_info->Nv, sizeof(PROROW));
}


void freeViewXYZInfo(struct ViewXYZInfo *view_xyz_info)
{
	free(view_xyz_info->ir_start);
	free(view_xyz_info->ir_num);
}

void createACol(struct ACol *A_col, int length)
{
	A_col->n_index =0;
	A_col->index = (int *)get_spc(length, sizeof(int));
	A_col->val = (ENTRY *)get_spc(length, sizeof(ENTRY));
	A_col->array_length = length;  /* sjk */
}

/* sjk */
void increaseAColLength(struct ACol *A_col)  /* increase array length by 10% */
{
        int *index,i,new_length;
        ENTRY *val;

        new_length=1.1*A_col->array_length;

        index = (int *)get_spc(new_length,sizeof(int));  /* create new array */
        for(i=0;i<A_col->array_length;i++)  /* copy old into new array */
                index[i] = A_col->index[i];
        free(A_col->index);  /* free old array */
        A_col->index = index;  /* point to new array */

        val = (ENTRY *)get_spc(new_length,sizeof(ENTRY));
        for(i=0;i<A_col->array_length;i++)
                val[i] = A_col->val[i];
        free(A_col->val);
        A_col->val= val;

        A_col->array_length= new_length;
}

void freeACol(struct ACol *A_col)
{
	free(A_col->index);
	free(A_col->val);
}


void compAColxyzOnFly(
	int jx,
	int jy,
	int jz,
	float x,
	float y,
	float z,
	struct GeomInfo *geom_info,
	struct SourceLocInfo *source_loc_info,
	struct ViewXYInfo *view_xy_info,
	struct ViewXYZInfo *view_xyz_info,
	struct ACol *A_col,struct ImgInfo *img_info)
{
	int iv, ic, ir, l, p, q;
	float del_r;
	float Cij;
	float Atmp;
	float inverse_num=1.0/255;

	A_col->n_index = 0;
	for (l = 0; l < view_xyz_info->iv_num; l++)	/* view loop */
	{
		iv = view_xyz_info->iv_start + l;
		float C_multiply = (sqrt((source_loc_info->xs[iv]-x)*(source_loc_info->xs[iv]-x) + (source_loc_info->ys[iv]-y)*(source_loc_info->ys[iv]-y) + (source_loc_info->zs[iv]-z)*(source_loc_info->zs[iv]-z)));
		float C_divide = (sqrt((source_loc_info->xs[iv]-x)*(source_loc_info->xs[iv]-x) + (source_loc_info->ys[iv]-y)*(source_loc_info->ys[iv]-y))*geom_info->Del_dr);
	
		for (p = 0; p < view_xy_info->ic_num[jx][jy][iv]; p++)	/* channel loop */
		{
			ic = view_xy_info->ic_start[jx][jy][iv] + p;
			for (q = 0; q < view_xyz_info->ir_num[iv]; q++)	/* row loop */
			{
				ir = view_xyz_info->ir_start[iv] + q;
				/* ATTENTION!! CHANGE SIGN HERE! ROW 0 IS CLOSEST */
				del_r = view_xy_info->Mag[jx][jy][iv]*view_xy_info->Mag_MaxPointer[jx*img_info->Ny+jy]*inverse_num*(z-source_loc_info->zs[iv]) + ir*geom_info->Del_dr - geom_info->half_detr;

				Cij = clip(0.0, (((img_info->Del_z*view_xy_info->Mag[jx][jy][iv]*view_xy_info->Mag_MaxPointer[jx*img_info->Ny+jy]*inverse_num+geom_info->Del_dr)/2.0)-fabs(del_r)), min(img_info->Del_z*view_xy_info->Mag[jx][jy][iv]*view_xy_info->Mag_MaxPointer[jx*img_info->Ny+jy]*inverse_num, geom_info->Del_dr));
				Cij *= C_multiply;
				Cij /= C_divide;
				Atmp = view_xy_info->B[jx][jy][iv][p]*view_xy_info->B_MaxPointer[jx*img_info->Ny+jy]*inverse_num*Cij;
				if (Atmp > EPSILON)	/* non-zero entry */
				{
					A_col->index[A_col->n_index] = iv*geom_info->Nc*geom_info->Nr + ic*geom_info->Nr + ir;
					A_col->val[A_col->n_index] = Atmp;
					A_col->n_index++;
					/* sjk: */
					if(A_col->n_index >= A_col->array_length)
					{
					  fprintf(stdout,"Increasing size of A column by 10%% (%d).\n",(int)1.1*A_col->array_length);
					  increaseAColLength(A_col);
					}
				}
			}
		}
	}
}






void metal_forward(ENTRY *e, ENTRY *X, char **recon_mask, struct GeomInfo *geom_info, struct ImgInfo *img_info,struct ViewXYInfo *view_xy_info, int myid,int total_nodes)
{
	int i, t;

	struct SourceLocInfo source_loc_info;

	/* initialize projection */
	memset(e, 0, (geom_info->Nr)*(geom_info->Nc)*(geom_info->Nv)*sizeof(ENTRY)); 	

	/* allocate precomputing structures */
	createSourceLocInfo(&source_loc_info, geom_info);
	compSourceLocInfo(&source_loc_info, geom_info,myid,total_nodes);
	

	/* for each voxel in the image */
	fprintf(stdout, "\nmetal forward projecting (parallelized version) ...\n");
	fflush(stdout);
	
	paraForwardProject(geom_info,img_info,&source_loc_info,X,e,recon_mask,view_xy_info);

	fprintf(stdout, "\nfinish metal forward projection!\n");
	fflush(stdout);


	freeSourceLocInfo(&source_loc_info);
}


void forwardProject(ENTRY *e, ENTRY *X, ENTRY *Y, char **recon_mask, struct GeomInfo *geom_info, struct ImgInfo *img_info,struct ViewXYInfo *view_xy_info, int myid,int total_nodes)
{
	int i, t;

	struct SourceLocInfo source_loc_info;

	/* initialize projection */
	memset(e, 0, (geom_info->Nr)*(geom_info->Nc)*(geom_info->Nv)*sizeof(ENTRY)); 	

	/* allocate precomputing structures */
	createSourceLocInfo(&source_loc_info, geom_info);
	compSourceLocInfo(&source_loc_info, geom_info,myid,total_nodes);
	

	/* for each voxel in the image */
	fprintf(stdout, "\nforward projecting (parallelized version) ...\n");
	fflush(stdout);
	
	paraForwardProject(geom_info,img_info,&source_loc_info,X,e,recon_mask,view_xy_info);

	for (i=0;i<(geom_info->Nr)*(geom_info->Nc)*(geom_info->Nv);i++){
		e[i] = Y[i]-e[i];
	}

	fprintf(stdout, "\nfinish forward projection!\n");



	freeSourceLocInfo(&source_loc_info);
}

void paraForwardProject(struct GeomInfo *geom_info,struct ImgInfo *img_info,struct SourceLocInfo *source_loc_info,ENTRY *X,ENTRY *e,char **recon_mask,	struct ViewXYInfo *view_xy_info)
{

	int offset1,offset2;  /* sjk */
	int jx, jy, jz,  Nyz;
	float x, y;


//	createViewXYInfo(view_xy_info, geom_info);
	createViewXYInfo_stored(view_xy_info, geom_info,img_info,recon_mask);

	


	Nyz = img_info->Ny*img_info->Nz;

	for (jx = 0; jx < img_info->Nx; jx++)
	{
		x = img_info->x0 + jx*img_info->Del_xy;
		offset1=jx*Nyz;
		
		for (jy = 0; jy < img_info->Ny; jy++)
		{
			if(recon_mask[jx][jy])
	
			{

				y = img_info->y0 + jy*img_info->Del_xy;
	
				
				compViewXYInfo(jx,jy,x, y, view_xy_info, geom_info, img_info, source_loc_info);

			

				offset2=offset1+jy*img_info->Nz;

				
				#pragma omp parallel for
				for (jz = 0; jz < img_info->Nz; jz++)
				{

					int indx=offset2+jz;  
					if(X[indx]>0) 
					{
						float z = img_info->z0 + jz*img_info->Del_z;
	
						struct ViewXYZInfo view_xyz_info;

						createViewXYZInfo(&view_xyz_info, geom_info);

						compViewXYZInfo(jx,jy,jz,z, &view_xyz_info, geom_info, img_info, source_loc_info, view_xy_info);

						struct ACol col_xyz;

						createACol(&col_xyz, COL_LEN);

						compAColxyzOnFly(jx,jy,jz,x, y, z, geom_info, source_loc_info, view_xy_info, &view_xyz_info, &col_xyz,img_info);
	
						for (int r = 0; r < col_xyz.n_index; r++)
						{
							#pragma omp atomic
							e[col_xyz.index[r]] += col_xyz.val[r]*X[indx];  
						}
						freeViewXYZInfo(&view_xyz_info);
					
						freeACol(&col_xyz);
					
					}  
				}
				
				
//				freeViewXYInfoB(view_xy_info, geom_info);
			}
		}
	}
	

//	freeViewXYInfo(view_xy_info, geom_info);

}


