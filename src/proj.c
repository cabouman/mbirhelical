#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sched.h>
#include <omp.h>
#include "allocate.h"
#include "data.h"
#include "prepro.h"
#include "proj.h"

void createSinogram(struct Sinogram *sinogram)
{
	sinogram->sino = (ENTRY *)get_spc((sinogram->geom_info.Nr)*(sinogram->geom_info.Nc)*(sinogram->geom_info.Nv), sizeof(ENTRY));
}

void freeSinogram(struct Sinogram *sinogram)
{
	free(sinogram->sino);
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
	geom_info->Del_zs = (geom_info->u)*(geom_info->Del_dr)*(geom_info->Del_beta)/(2.0*PI);
	geom_info->detc = (geom_info->Nc)*(geom_info->Del_alphac);
	geom_info->detr = (geom_info->Nr)*(geom_info->Del_dr);
	geom_info->half_detr = (geom_info->Nr-1.0)*(geom_info->Del_dr)/2.0;
	geom_info->cone_zbuffer= geom_info->detr*(geom_info->r_si+geom_info->fov/2.0)/(2.0*geom_info->r_sd);
}

void fillImgInfo(struct ImgInfo *img_info)	/* fill in the intermediate variables */
{
	img_info->x0 = img_info->xc - (img_info->Del_xy)*(img_info->Nx-1)/2.0;
	img_info->y0 = img_info->yc - (img_info->Del_xy)*(img_info->Ny-1)/2.0;
	img_info->z0 = img_info->zc - (img_info->Del_z)*(img_info->Nz-1)/2.0;
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
	struct GeomInfo *geom_info)
{
	int iv;		/* view index */
	for (iv = 0; iv < geom_info->Nv; iv++)
	{
		source_loc_info->beta[iv] = geom_info->beta0 + iv*(geom_info->Del_beta);
		source_loc_info->xs[iv] = (geom_info->r_si)*cos(source_loc_info->beta[iv]);
		source_loc_info->ys[iv] = (geom_info->r_si)*sin(source_loc_info->beta[iv]);
		source_loc_info->zs[iv] = iv*(geom_info->Del_zs);		/* assume z start from 0 */
	}
}

void freeSourceLocInfo(struct SourceLocInfo *source_loc_info)
{
	free(source_loc_info->beta);
	free(source_loc_info->xs);
	free(source_loc_info->ys);
	free(source_loc_info->zs);
}

void createViewXYInfo(
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info)
{
	view_xy_info->ic_start = (CHANNEL *)get_spc(geom_info->Nv, sizeof(CHANNEL));
	view_xy_info->ic_num = (PROCHANNEL *)get_spc(geom_info->Nv, sizeof(PROCHANNEL));
	view_xy_info->Mag = (ENTRY *)get_spc(geom_info->Nv, sizeof(ENTRY));
	view_xy_info->Wr = (ENTRY *)get_spc(geom_info->Nv, sizeof(ENTRY));

	view_xy_info->B = (ENTRY **)get_spc(geom_info->Nv, sizeof(ENTRY *));
}

void compViewXYInfo(
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
		view_xy_info->Mag[iv] = (geom_info->r_sd)/r_sv;
		view_xy_info->Wr[iv] = (img_info->Del_z)*(view_xy_info->Mag[iv]);

		Wc = (img_info->Del_xy)*costh/r_sv;


		alpha_min = alphaj_td - geom_info->alphac0 - (Wc - geom_info->Del_alphac)/2.0;
		alpha_max = alphaj_td - geom_info->alphac0 + (Wc + geom_info->Del_alphac)/2.0;

		if (alpha_max < 0 || alpha_min > geom_info->detc)
		{
			view_xy_info->ic_num[iv] = 0;
		}
		else
		{
			view_xy_info->ic_start[iv] = (CHANNEL)max((CHANNEL)floor(alpha_min/(geom_info->Del_alphac)), 0);
			ic_end = (CHANNEL)min((CHANNEL)floor(alpha_max/(geom_info->Del_alphac)), (CHANNEL)(geom_info->Nc-1));
			view_xy_info->ic_num[iv] = ((PROCHANNEL)(ic_end - view_xy_info->ic_start[iv] + 1));
		}
		

		
		view_xy_info->B[iv] = (ENTRY *)get_spc((int)(view_xy_info->ic_num[iv]), sizeof(ENTRY));

		
		for (p = 0; p < view_xy_info->ic_num[iv]; p++)
		{
			ic = (int)(view_xy_info->ic_start[iv] + p);

			del_c = adjust(alphaj - (ic*(geom_info->Del_alphac) + geom_info->alphac0));

			Bij = clip(0.0, ((Wc+(geom_info->Del_alphac))/2.0)-fabs(del_c), min(Wc, (geom_info->Del_alphac)));
			Bij *= ((img_info->Del_xy)/((geom_info->Del_alphac)*costh));

			view_xy_info->B[iv][p] = Bij;
		}
	

	}
}

void freeViewXYInfoB(
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info)
{
	int iv;		/* view index */
	for (iv = 0; iv < geom_info->Nv; iv++)
	{
		free(view_xy_info->B[iv]);
	}
}

void freeViewXYInfo(
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info)
{
	free(view_xy_info->ic_start);
	free(view_xy_info->ic_num);
	free(view_xy_info->Mag);
	free(view_xy_info->Wr);
	free(view_xy_info->B);
}

void createViewXYZInfo(
	struct ViewXYZInfo *view_xyz_info,
	struct GeomInfo *geom_info)
{
	view_xyz_info->ir_start = (ROW *)get_spc(geom_info->Nv, sizeof(ROW));
	view_xyz_info->ir_num = (PROROW *)get_spc(geom_info->Nv, sizeof(PROROW));
}

void compViewXYZInfo(
	float z,
	struct ViewXYZInfo *view_xyz_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info,
	struct SourceLocInfo *source_loc_info,
	struct ViewXYInfo *view_xy_info)
{
	int iv;		/* view index */
	float d;
	float d_min;
	float d_max;
	ROW ir_end;
	int iv_end;  /* sjk */


  if(0)  /* sjk: this block is replaced with the lines following */
  {
	/* calculate d */
	d = geom_info->Nr*geom_info->Del_dr*(geom_info->r_si+geom_info->fov/2.0)/(2.0*geom_info->r_sd);
	if (z < source_loc_info->zs[0]-d || z > source_loc_info->zs[geom_info->Nv-1]+d)
	{
		view_xyz_info->iv_num = 0;
	}
	else
	{
		for (iv = 0; iv < geom_info->Nv; iv++)
		{
			if (z >= source_loc_info->zs[iv]-d && z <= source_loc_info->zs[iv]+d)
			{
				view_xyz_info->iv_start = iv;
				break;
			}
		}
		view_xyz_info->iv_num = 0;
		for (iv = 0; iv < geom_info->Nv; iv++)
		{
			if (z >= source_loc_info->zs[iv]-d && z <= source_loc_info->zs[iv]+d)
			{
				view_xyz_info->iv_num++;
			}
		}
		/*for (; z >= source_loc_info->zs[iv]-d && z <= source_loc_info->zs[iv]+d; iv++)
		{
			view_xyz_info->iv_num++;
		}*/
	}
  }


	/* sjk: this block replaces the above, finding the iv range in closed form */
        view_xyz_info->iv_start=max(0,(int)floor((z-img_info->Del_z/2.0-geom_info->cone_zbuffer)/geom_info->Del_zs));
        iv_end =   min(geom_info->Nv-1,(int)ceil((z+img_info->Del_z/2.0+geom_info->cone_zbuffer)/geom_info->Del_zs));
        if((iv_end<0) || (view_xyz_info->iv_start > (geom_info->Nv-1)))
                view_xyz_info->iv_num=0;
        else
                view_xyz_info->iv_num= iv_end - view_xyz_info->iv_start + 1;

	/* sjk: moved this block down so that we don't have to go through all views */
	for (iv = view_xyz_info->iv_start; iv <= iv_end; iv++)   /* sjk */
	{
		d_min = (geom_info->half_detr+geom_info->Del_dr/2.0-view_xy_info->Mag[iv]*(z-source_loc_info->zs[iv]+img_info->Del_z/2.0));
		d_max = (geom_info->half_detr+geom_info->Del_dr/2.0-view_xy_info->Mag[iv]*(z-source_loc_info->zs[iv]-img_info->Del_z/2.0));

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

void freeViewXYZInfo(struct ViewXYZInfo *view_xyz_info)
{
	free(view_xyz_info->ir_start);
	free(view_xyz_info->ir_num);
}

void createACol(struct ACol *A_col, int length)
{
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

void compAColxyz(
	int jjx,
	int jjy,
	int jjz,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info,
	struct ACol *col_xyz)
{
	static int first = 1;
	static struct ACol A_col;
	static struct ACol *A_col_pt;
	static struct ACol ***A_array;		/* store the whole A matrix, each column is A[jx][jy][jz] */

	int jx, jy, jz;			/* voxel indicies */
	int iv, ic, ir, l, p, q;	/* detector indicies */
	int r;
	float x, y, z;			/* voxel coordinate */
	float del_r;
	float Atmp, Cij;
	struct SourceLocInfo source_loc_info;
	struct ViewXYInfo view_xy_info;
	struct ViewXYZInfo view_xyz_info;

	if (first == 1)			/* if the function is first called, precompute the whole A matrix and store it */
	{
		first = 0;

		A_array = (struct ACol ***)multialloc(sizeof(struct ACol), 3, img_info->Nx, img_info->Ny, img_info->Nz);/* no allocation for arrays */

		createACol(&A_col, COL_LEN);		/* TODO, COL_LEN hard-coded */

		/* allocate precomputing structure */
		createSourceLocInfo(&source_loc_info, geom_info);
		compSourceLocInfo(&source_loc_info, geom_info);		/* populate SourceLocInfo structure */
		createViewXYInfo(&view_xy_info, geom_info);
		createViewXYZInfo(&view_xyz_info, geom_info);

		fprintf(stdout, "\ninitializing A matrix ...\n");

		for (jx = 0 ; jx < img_info->Nx; jx++)
		{
			x = img_info->x0 + jx*img_info->Del_xy;
			for (jy = 0; jy < img_info->Ny; jy++)
			{
				y = img_info->y0 + jy*img_info->Del_xy;

				compViewXYInfo(x, y, &view_xy_info, geom_info, img_info, &source_loc_info);	/* populate ViewXYInfo structure */

				for (jz = 0; jz < img_info->Nz; jz++)
				{
					z = img_info->z0 + jz*img_info->Del_z;

					compViewXYZInfo(z, &view_xyz_info, geom_info, img_info, &source_loc_info, &view_xy_info);	/* poplutate ViewXYZInfo structure */

					A_col.n_index = 0;
					for (l = 0; l < view_xyz_info.iv_num; l++)	/* view loop */
					{
						iv = view_xyz_info.iv_start + l;
						for (p = 0; p < view_xy_info.ic_num[iv]; p++)	/* channel loop */
						{
							ic = view_xy_info.ic_start[iv] + p;
							for (q = 0; q < view_xyz_info.ir_num[iv]; q++)	/* row loop */
							{
								ir = view_xyz_info.ir_start[iv] + q;
								/* ATTENTION! CHANGE SIGN HERE ROW 0 IS CLOSEST !! */
								del_r = view_xy_info.Mag[iv]*(z-source_loc_info.zs[iv]) + ir*geom_info->Del_dr - geom_info->half_detr;
								Cij = clip(0.0, ((view_xy_info.Wr[iv]+geom_info->Del_dr)/2.0)-fabs(del_r), min(view_xy_info.Wr[iv], geom_info->Del_dr));  /* sjk: fixed typo */
								Cij *= (sqrt((source_loc_info.xs[iv]-x)*(source_loc_info.xs[iv]-x)+(source_loc_info.ys[iv]-y)*(source_loc_info.ys[iv]-y)+(source_loc_info.zs[iv]-z)*(source_loc_info.zs[iv]-z)));
								Cij /= (sqrt((source_loc_info.xs[iv]-x)*(source_loc_info.xs[iv]-x) +(source_loc_info.ys[iv]-y)*(source_loc_info.ys[iv]-y))*geom_info->Del_dr);
								Atmp = view_xy_info.B[iv][p]*Cij;
								if (Atmp > EPSILON)	/* if nonzero entry */
								{
									A_col.index[A_col.n_index] = iv*geom_info->Nc*geom_info->Nr + ic*geom_info->Nr + ir;
									A_col.val[A_col.n_index] = Atmp;
									A_col.n_index++;
								}
							}
						}
					}

					/* here we finish computing one column of A for a specific (x,y,z) voxel */
					/* store it in A_array */

					A_array[jx][jy][jz].index = (int *)get_spc(A_col.n_index, sizeof(int));
					A_array[jx][jy][jz].val = (ENTRY *)get_spc(A_col.n_index, sizeof(ENTRY));
					A_array[jx][jy][jz].n_index = A_col.n_index;
					for (r = 0; r < A_col.n_index; r++)
					{
						A_array[jx][jy][jz].index[r] = (int)A_col.index[r];
						A_array[jx][jy][jz].val[r] = (ENTRY)A_col.val[r];
					}
				}
				freeViewXYInfoB(&view_xy_info, geom_info);
			}
		}
		freeViewXYZInfo(&view_xyz_info);
		freeViewXYInfo(&view_xy_info, geom_info);
		freeSourceLocInfo(&source_loc_info);

		fprintf(stdout, "finish computing A matrix!\n");
	}

	A_col_pt = &A_array[jjx][jjy][jjz];
	col_xyz->n_index = A_col_pt->n_index;
	for (r = 0; r < A_col_pt->n_index; r++)
	{
		col_xyz->index[r] = A_col_pt->index[r];
		col_xyz->val[r] = A_col_pt->val[r];
	}
}

void compAColxyzOnFly(
	float x,
	float y,
	float z,
	struct GeomInfo *geom_info,
	struct SourceLocInfo *source_loc_info,
	struct ViewXYInfo *view_xy_info,
	struct ViewXYZInfo *view_xyz_info,
	struct ACol *A_col)
{
	int iv, ic, ir, l, p, q;
	float del_r;
	float Cij;
	float Atmp;

	A_col->n_index = 0;
	for (l = 0; l < view_xyz_info->iv_num; l++)	/* view loop */
	{
		iv = view_xyz_info->iv_start + l;
		for (p = 0; p < view_xy_info->ic_num[iv]; p++)	/* channel loop */
		{
			ic = view_xy_info->ic_start[iv] + p;
			for (q = 0; q < view_xyz_info->ir_num[iv]; q++)	/* row loop */
			{
				ir = view_xyz_info->ir_start[iv] + q;
				/* ATTENTION!! CHANGE SIGN HERE! ROW 0 IS CLOSEST */
				del_r = view_xy_info->Mag[iv]*(z-source_loc_info->zs[iv]) + ir*geom_info->Del_dr - geom_info->half_detr;
				Cij = clip(0.0, (((view_xy_info->Wr[iv]+geom_info->Del_dr)/2.0)-fabs(del_r)), min(view_xy_info->Wr[iv], geom_info->Del_dr));
				Cij *= (sqrt((source_loc_info->xs[iv]-x)*(source_loc_info->xs[iv]-x) + (source_loc_info->ys[iv]-y)*(source_loc_info->ys[iv]-y) + (source_loc_info->zs[iv]-z)*(source_loc_info->zs[iv]-z)));
				Cij /= (sqrt((source_loc_info->xs[iv]-x)*(source_loc_info->xs[iv]-x) + (source_loc_info->ys[iv]-y)*(source_loc_info->ys[iv]-y))*geom_info->Del_dr);
				Atmp = view_xy_info->B[iv][p]*Cij;
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

void forwardProject(ENTRY *AX, ENTRY *X, unsigned short *AX_mask, char **recon_mask, struct GeomInfo *geom_info, struct ImgInfo *img_info)
// Computes the error sinogram update e = Y - AX
{
	int i, t;

	struct SourceLocInfo source_loc_info;

	/* initialize projection */
	for (i = 0; i < (geom_info->Nr)*(geom_info->Nc)*(geom_info->Nv); i++)
	{
		AX[i] = 0.0;
	}
	if(AX_mask != NULL)
	{
		for (i = 0; i < (geom_info->Nr)*(geom_info->Nc)*(geom_info->Nv); i++)
			AX_mask[i] = 0;
	}

	/* allocate precomputing structures */
	createSourceLocInfo(&source_loc_info, geom_info);
	compSourceLocInfo(&source_loc_info, geom_info);
	

	/* for each voxel in the image */
	fprintf(stdout, "\nforward projecting (parallelized version) ...\n");

	
	#pragma omp parallel
	{
		paraForwardProject(geom_info,img_info,&source_loc_info,X,AX,AX_mask,recon_mask);
	}

	fprintf(stdout, "\nfinish forward projection!\n");



	freeSourceLocInfo(&source_loc_info);
}

void paraForwardProject(struct GeomInfo *geom_info,struct ImgInfo *img_info,struct SourceLocInfo *source_loc_info,ENTRY *X,ENTRY *AX,unsigned short *AX_mask,char **recon_mask)
// Computes the forward projection of X and returns in e
// A matrix is recomputed column-by-column every call
{
	int tid = omp_get_thread_num();

	int indx,offset1,offset2;  /* sjk */
	int jx, jy, jz, jzmax, Nyz, r;
	float x, y, z;
	struct ACol col_xyz;
	struct ViewXYInfo view_xy_info;
	struct ViewXYZInfo view_xyz_info;


	createACol(&col_xyz, COL_LEN);		/* TODO COL_LEN hard-coded */
	createViewXYInfo(&view_xy_info, geom_info);
	createViewXYZInfo(&view_xyz_info, geom_info);


	fprintf(stdout,"tid is %d omp_get_num_threads() %d\n",tid,omp_get_num_threads());

	jzmax = (tid+1)*img_info->Nz/omp_get_num_threads();
	if (tid == (omp_get_num_threads()-1))
	{
		jzmax = img_info->Nz;
	}
	Nyz = img_info->Ny*img_info->Nz;

	for (jx = 0; jx < img_info->Nx; jx++)
	{
		x = img_info->x0 + jx*img_info->Del_xy;
		offset1=jx*Nyz;
		
		for (jy = 0; jy < img_info->Ny; jy++)
		{
			if(recon_mask[jx][jy])
		   /* sjk: everything outside mask has been set to 0, and never updated */
	
			{
			
				y = img_info->y0 + jy*img_info->Del_xy;
	
				
				compViewXYInfo(x, y, &view_xy_info, geom_info, img_info, source_loc_info);
				

				
				offset2=offset1+jy*img_info->Nz; 
				for (jz = tid*img_info->Nz/omp_get_num_threads(); jz < jzmax; jz++)
				{
					indx=offset2+jz;  
					if(X[indx]>0) 
					{
						z = img_info->z0 + jz*img_info->Del_z;
		
						compViewXYZInfo(z, &view_xyz_info, geom_info, img_info, source_loc_info, &view_xy_info);
							
						compAColxyzOnFly(x, y, z, geom_info, source_loc_info, &view_xy_info, &view_xyz_info, &col_xyz);
	
						#pragma omp critical
						{
							for (r = 0; r < col_xyz.n_index; r++)
							{
								AX[col_xyz.index[r]] += col_xyz.val[r]*X[indx];  
								if(AX_mask != NULL)
								{
									if(X[indx]>hu2miu(3000,MIU_AIR,MIU_WATER)) 
										AX_mask[col_xyz.index[r]] += 1; 
								}

							}
						}
						
						
					}  
				}
				
				
				freeViewXYInfoB(&view_xy_info, geom_info);
				
			}
		}
	}
	

	freeViewXYZInfo(&view_xyz_info);
	freeViewXYInfo(&view_xy_info, geom_info);
	freeACol(&col_xyz);

}

void serialForwardProject(ENTRY *AX, ENTRY *X, struct GeomInfo *geom_info, struct ImgInfo *img_info)
{
	int i, jx, jy, jz, r, Nyz;
	float x, y, z;

	struct SourceLocInfo source_loc_info;
	struct ViewXYInfo view_xy_info;
	struct ViewXYZInfo view_xyz_info;
	struct ACol col_xyz;

	/* initialize projection */
	for (i = 0; i < (geom_info->Nr)*(geom_info->Nc)*(geom_info->Nv); i++)
	{
		AX[i] = 0.0;
	}

	/* allocate precomputing structures */
	createACol(&col_xyz, COL_LEN);		/* TODO COL_LEN hard-coded */
	createSourceLocInfo(&source_loc_info, geom_info);
	compSourceLocInfo(&source_loc_info, geom_info);
	createViewXYInfo(&view_xy_info, geom_info);
	createViewXYZInfo(&view_xyz_info, geom_info);

	Nyz = img_info->Ny*img_info->Nz;

	/* for each voxel in the image */
	fprintf(stdout, "\nforward projecting (serial version) ...\n");
	for (jx = 0; jx < img_info->Nx; jx++)
	{
		x = img_info->x0 + jx*img_info->Del_xy;
		for (jy = 0; jy < img_info->Ny; jy++)
		{
			y = img_info->y0 + jy*img_info->Del_xy;

			compViewXYInfo(x, y, &view_xy_info, geom_info, img_info, &source_loc_info);

			for (jz = 0; jz < img_info->Nz; jz++)
			{
				z = img_info->z0 + jz*img_info->Del_z;

				compViewXYZInfo(z, &view_xyz_info, geom_info, img_info, &source_loc_info, &view_xy_info);

				compAColxyzOnFly(x, y, z, geom_info, &source_loc_info, &view_xy_info, &view_xyz_info, &col_xyz);
				for (r = 0; r < col_xyz.n_index; r++)
				{
					AX[col_xyz.index[r]] += col_xyz.val[r]*X[jx*Nyz+jy*img_info->Nz+jz];
				}
			}
			freeViewXYInfoB(&view_xy_info, geom_info);
		}
	}
	fprintf(stdout, "\nfinish forward projection!\n");

	freeViewXYZInfo(&view_xyz_info);
	freeViewXYInfo(&view_xy_info, geom_info);
	freeSourceLocInfo(&source_loc_info);
	freeACol(&col_xyz);
}




void backProject(ENTRY *AX, ENTRY *X, struct GeomInfo *geom_info, struct ImgInfo *img_info)
{
	int i, t, len;

	struct SourceLocInfo source_loc_info;
	//pthread_t thread[NUM_CORE];
	//struct paraForwardProjectData thread_data[NUM_CORE];

	/* initialize projection */
	len=img_info->Nx * img_info->Ny * img_info->Nz;
	for (i = 0; i < len; i++)
		X[i] = 0.0;

	/* allocate precomputing structures */
	createSourceLocInfo(&source_loc_info, geom_info);
	compSourceLocInfo(&source_loc_info, geom_info);

	/* for each voxel in the image */
	fprintf(stdout, "\nback projecting (parallelized version) ...\n");
	/*
	for (t = 0; t < NUM_CORE; t++)
	{
		thread_data[t].tid = t;
		thread_data[t].geom_info = geom_info;
		thread_data[t].img_info = img_info;
		thread_data[t].source_loc_info = &source_loc_info;
		thread_data[t].X = X;
		thread_data[t].AX = AX;

		pthread_create(&thread[t], NULL, paraBackProject, (void *)&thread_data[t]);
	}
	for (t = 0; t < NUM_CORE; t++)
	{
		pthread_join(thread[t], NULL);
	}

	*/
	paraBackProject(geom_info,img_info,&source_loc_info,X,AX);

	freeSourceLocInfo(&source_loc_info);
}


void *paraBackProject(struct GeomInfo *geom_info,struct ImgInfo *img_info,struct SourceLocInfo *source_loc_info,ENTRY *X,ENTRY *AX)
{
	int tid = omp_get_thread_num();
	//struct GeomInfo *geom_info;
	//struct ImgInfo *img_info;
	//struct SourceLocInfo *source_loc_info;
	//ENTRY *X;
	//ENTRY *AX;
	//struct paraForwardProjectData *data;

	int indx,offset1,offset2;  /* sjk */
	int jx, jy, jz, jzmax, Nyz, r;
	float x, y, z;
	struct ACol col_xyz;
	struct ViewXYInfo view_xy_info;
	struct ViewXYZInfo view_xyz_info;
	float sum;

	//data = (struct paraForwardProjectData *)input;
	//tid = data->tid;
	//geom_info = data->geom_info;
	//img_info = data->img_info;
	//source_loc_info = data->source_loc_info;
	//X = data->X;
	//AX = data->AX;

	createACol(&col_xyz, COL_LEN);		/* TODO COL_LEN hard-coded */
	createViewXYInfo(&view_xy_info, geom_info);
	createViewXYZInfo(&view_xyz_info, geom_info);

	jzmax = (tid+1)*img_info->Nz/omp_get_num_threads();
	if (tid == (omp_get_num_threads()-1))
	{
		jzmax = img_info->Nz;
	}
	Nyz = img_info->Ny*img_info->Nz;

	for (jx = 0; jx < img_info->Nx; jx++)
	{
		x = img_info->x0 + jx*img_info->Del_xy;
		offset1=jx*Nyz;  /* sjk */
		for (jy = 0; jy < img_info->Ny; jy++)
		{
			y = img_info->y0 + jy*img_info->Del_xy;

			compViewXYInfo(x, y, &view_xy_info, geom_info, img_info, source_loc_info);

			offset2=offset1+jy*img_info->Nz;  /* sjk */
			for (jz = tid*img_info->Nz/omp_get_num_threads(); jz < jzmax; jz++)
			{
				indx=offset2+jz;  /* sjk */
				z = img_info->z0 + jz*img_info->Del_z;

				compViewXYZInfo(z, &view_xyz_info, geom_info, img_info, source_loc_info, &view_xy_info);

				compAColxyzOnFly(x, y, z, geom_info, source_loc_info, &view_xy_info, &view_xyz_info, &col_xyz);
				sum=0;
				for (r = 0; r < col_xyz.n_index; r++)
				{
					sum += AX[col_xyz.index[r]] * col_xyz.val[r]; 
					/*AX[col_xyz.index[r]] += col_xyz.val[r]*X[indx]; */ 
				}
				X[indx]=sum/col_xyz.n_index;  /* divide by number of non-zero terms */
			}
			freeViewXYInfoB(&view_xy_info, geom_info);
		}
	}

	freeViewXYZInfo(&view_xyz_info);
	freeViewXYInfo(&view_xy_info, geom_info);
	freeACol(&col_xyz);

	return 0;
}

