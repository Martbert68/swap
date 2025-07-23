#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#include "martin.h"
#define BRIDGE 50
#define DOWN 137
#define UP 250

/* here are our X variables */
Display *dis;

int screen,len;
double scale;
double mass[6];
double damp;
int fret[6];
unsigned char *x_buffer;
Window win;
GC gc;
XImage *x_image;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();
int f_col,f_row,f_off,dot_x,dot_y;
static double pp; 



void disp (double *,int,int);

void usage ()
{
	printf("usage: font filename threshold [20-40 ish] star,framestcode [65A 97a]\n");
	exit (1);
}

int qlott(unsigned char *image, int xp, int yp, int r, int g, int b, int xs, int ys)
{
	int x,y;
	if (xs<1){ xs=1;}
	if (ys<1){ ys=1;}
	for (x=0;x<xs;x++)
	{
		for (y=0;y<ys;y++)
		{
			plott(image,xp+x,yp+y,r,g,b);
		}
	}
}	



double score (double *image, int x, int y)
{
	int ul,um,ur;
	int ml,mm,mr;
	int dl,dm,dr;

	double totalr,totalg,totalb;
	double avgr,avgg,avgb,sr,sg,sb,sd;

	ul=((x-1)*3)+((y-1)*3*X_SIZE);
	um=(x*3)+((y-1)*3*X_SIZE);
	ur=((x+1)*3)+((y-1)*3*X_SIZE);

	ml=((x-1)*3)+(y*3*X_SIZE);
	mm=(x*3)+(y*3*X_SIZE);
	mr=((x+1)*3)+(y*3*X_SIZE);

	dl=((x-1)*3)+((y+1)*3*X_SIZE);
	dm=(x*3)+((y+1)*3*X_SIZE);
	dr=((x+1)*3)+((y+1)*3*X_SIZE);


	totalr=image[ul]+image[um]+image[ur]+image[ml]+image[mm]+image[mr]+image[dl]+image[dm]+image[dr] ;
	totalg=image[ul+1]+image[um+1]+image[ur+1]+image[ml+1]+image[mm+1]+image[mr+1]+image[dl+1]+image[dm+1]+image[dr+1] ;
	totalb=image[ul+2]+image[um+2]+image[ur+2]+image[ml+2]+image[mm+2]+image[mr+2]+image[dl+2]+image[dm+2]+image[dr+2] ;

	avgr=totalr/9;	
	avgg=totalg/9;	
	avgb=totalb/9;	

	if (avgr==0){avgr=1;}
	if (avgg==0){avgg=1;}
	if (avgb==0){avgb=1;} 

	sr=((avgr-image[ul])*(avgr-image[ul]))+((avgr-image[ur])*(avgr-image[ur]))+((avgr-image[um])*(avgr-image[um]))+
	   ((avgr-image[ml])*(avgr-image[ml]))+((avgr-image[mr])*(avgr-image[mr]))+((avgr-image[mm])*(avgr-image[mm]))+
	   ((avgr-image[dl])*(avgr-image[dl]))+((avgr-image[dr])*(avgr-image[dr]))+((avgr-image[dm])*(avgr-image[dm]));

	sg=((avgg-image[ul+1])*(avgg-image[ul+1]))+((avgg-image[ur+1])*(avgg-image[ur+1]))+((avgg-image[um+1])*(avgg-image[um+1]))+
	   ((avgg-image[ml+1])*(avgg-image[ml+1]))+((avgg-image[mr+1])*(avgg-image[mr+1]))+((avgg-image[mm+1])*(avgg-image[mm+1]))+
	   ((avgg-image[dl+1])*(avgg-image[dl+1]))+((avgg-image[dr+1])*(avgg-image[dr+1]))+((avgg-image[dm+1])*(avgg-image[dm+1]));


	sb=((avgb-image[ul+2])*(avgb-image[ul+2]))+((avgb-image[ur+2])*(avgb-image[ur+2]))+((avgb-image[um+2])*(avgb-image[um+2]))+
	   ((avgb-image[ml+2])*(avgb-image[ml+2]))+((avgb-image[mr+2])*(avgb-image[mr+2]))+((avgb-image[mm+2])*(avgb-image[mm+2]))+
	   ((avgb-image[dl+2])*(avgb-image[dl+2]))+((avgb-image[dr+2])*(avgb-image[dr+2]))+((avgb-image[dm+2])*(avgb-image[dm+2]));

	sd=sqrt(sr)+sqrt(sg)+sqrt(sb);
	//sd=sr+sg+sb;

	return sd;

}



void swap (double *image,int x, int y, int dxa, int dya, int dxb, int dyb)
{
	double sr,sg,sb;
	double dr,dg,db;

	sr=image[((x)*3)+(3*(y)*X_SIZE)];
	sg=image[((x)*3)+(3*(y)*X_SIZE)+1];
	sb=image[((x)*3)+(3*(y)*X_SIZE)+2];

	dr=image[((dxb)*3)+(3*(dyb)*X_SIZE)];
	dg=image[((dxb)*3)+(3*(dyb)*X_SIZE)+1];
	db=image[((dxb)*3)+(3*(dyb)*X_SIZE)+2];

	image[((x)*3)+(3*(y)*X_SIZE)]=dr;
	image[((x)*3)+(3*(y)*X_SIZE)+1]=dg;
	image[((x)*3)+(3*(y)*X_SIZE)+2]=db;

	image[((dxb)*3)+(3*(dyb)*X_SIZE)]=sr;
	image[((dxb)*3)+(3*(dyb)*X_SIZE)+1]=sg;
	image[((dxb)*3)+(3*(dyb)*X_SIZE)+2]=sb;

}


void lineq (double *image, int sx, int sy, int ex, int ey)
{
	int x,y,along,len,dx,dy;
	int rs,re,rd,gs,ge,gd,bs,be,bd;

	rs=image[(sx*3)+(sy*3*X_SIZE)+0]; re=image[(ex*3)+(ey*3*X_SIZE)+0]; rd=re-rs;
	gs=image[(sx*3)+(sy*3*X_SIZE)+1]; ge=image[(ex*3)+(ey*3*X_SIZE)+1]; gd=ge-gs;
	bs=image[(sx*3)+(sy*3*X_SIZE)+2]; be=image[(ex*3)+(ey*3*X_SIZE)+2]; bd=be-bs;

	dx=ex-sx; dy=ey-sy;

	len=sqrt((dx*dx)+(dy*dy));

	for (along=0;along<len;along+=13)
	{
		x=(along*dx/len)+sx;
		y=(along*dy/len)+sy;


		image[(x*3)+(y*X_SIZE*3)]=rs+(along*rd/len);
		image[(x*3)+(y*X_SIZE*3)+1]=gs+(along*gd/len);
		image[(x*3)+(y*X_SIZE*3)+2]=bs+(along*bd/len);

	}
}




int main(int argc,char *argv[])
{
	double *image;
	int i,j,loop,frame;
	unsigned char *image1;
	unsigned char *imagel;

        image=(double *)malloc(sizeof (double)*X_SIZE*Y_SIZE*3); // disp buffer
        image1=(char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        imagel=(char *)malloc(sizeof (char)*X_SIZE*Y_SIZE); //locker 


	char junk[30];


	init_x();

	frame=0;

	for (i=0;i<X_SIZE*Y_SIZE;i++){ imagel[i]=0;}

	for (i=0;i<3*X_SIZE*Y_SIZE;i+=3){ 

		double intensity,red,green,blue,tot,r,g,b;

		intensity=74+(double)(rand()%712000)/1000;

		
		red=(double)(rand()%10000)/10000;
		green=(double)(rand()%10000)/10000;
		blue=(double)(rand()%10000)/10000;
		tot=red+green+blue;
		if (tot==0){ tot=1;}
		
		r=red*intensity/tot;
		g=green*intensity/tot;
		b=blue*intensity/tot;

		if (r>255){r=255;}
		if (g>255){g=255;}
		if (b>255){b=255;}

		image[i]=r;
		image[i+1]=g;
		image[i+2]=b;

	}
	/*
	for (i=0;i<3*X_SIZE*Y_SIZE;i++){ image[i]=rand()%256;}
	for (i=0;i<X_SIZE;i++)
	{ 
		for (j=0;j<Y_SIZE;j++)
		{
			image[(i*3)+(j*3*X_SIZE)+0]=128+(127*sin((double)(((i*i))+(j*j))/26050));
			image[(i*3)+(j*3*X_SIZE)+1]=128+(127*sin((double)(((i*i))+(j*j))/33070));
			image[(i*3)+(j*3*X_SIZE)+2]=128+(127*sin((double)(((i*i))+(j*j))/57900));
		}
	} 

	for (i=0;i<3*X_SIZE*Y_SIZE;i++){ image[i]=20+rand()%200;} */

	load_image (image1,"sc.jpg",3);
	load_image (image1,"sc.jpg",0);
	for (i=0;i<3*X_SIZE*Y_SIZE;i++){ image[i]=image1[i];} 
	//for (i=0;i<3*X_SIZE*Y_SIZE;i++){ image[i]=(double)(rand()%255001)/1000;} 
	int xx,yy;
	for (xx=0;xx<X_SIZE;xx++)
	{
		for (yy=0;yy<Y_SIZE;yy++)
		{
			if (xx%90<30){  
				image[(xx*3)+(yy*3*X_SIZE)]=200+rand()%5;
				image[(xx*3)+(yy*3*X_SIZE)+1]=rand()%5;
				image[(xx*3)+(yy*3*X_SIZE)+2]=rand()%5; }
			else if (xx%93<61){
				image[(xx*3)+(yy*3*X_SIZE)]=rand()%5;
				image[(xx*3)+(yy*3*X_SIZE)+1]=200+rand()%5;
				image[(xx*3)+(yy*3*X_SIZE)+2]=rand()%5; }
			else {
				image[(xx*3)+(yy*3*X_SIZE)]=0;
				image[(xx*3)+(yy*3*X_SIZE)+1]=0;
				image[(xx*3)+(yy*3*X_SIZE)+2]=200+rand()%5; }

			if (yy%93<35){  image[(xx*3)+(yy*3*X_SIZE)+2]=200+rand()%5;}
				
			else if (yy%91<67){
				image[(xx*3)+(yy*3*X_SIZE)+1]=200+rand()%5;}
			else {
				image[(xx*3)+(yy*3*X_SIZE)]=200+rand()%5;}

		}
	} 


	double after;
	double before;
	int x,y,xa,ya;
	long int f;
	f=0;


	x=4+rand()%(X_SIZE-8);y=4+rand()%(Y_SIZE-8);
	x=4+rand()%(X_SIZE-8);y=4+rand()%(Y_SIZE-8);

	x=1;y=1;


	while (1)
	{
		

		//for (x=4;x<X_SIZE-4;x++){
		//for (y=4;y<Y_SIZE-4;y++){
		x=1+rand()%(X_SIZE-2);y=1+rand()%(Y_SIZE-2);
		xa=1+rand()%(X_SIZE-2);ya=1+rand()%(Y_SIZE-2);
		xa=1+rand()%(X_SIZE-2);ya=1+rand()%(Y_SIZE-2);

		xa=x-80+rand()%161;
		ya=y-80+rand()%161;

		if (xa<0 || xa >= X_SIZE){ continue;}
		if (ya<0 || ya >= Y_SIZE){ continue;} 

		if (imagel[x+(y*X_SIZE)]>0){imagel[x+(y*X_SIZE)]-- ;continue;}


		double dx,dy;

		dx=x-xa;dy=y-ya;

		if ((dx*dx)+(dy*dy)>4900){ continue;}
		

		/*double dx,dy;

		dx=(x-xa);
		dy=(y-ya);
		if ((dx*dx)+(dy*dy)>400){ continue;}*/


		double aa,bb;

		aa=score(image,x,y); 
		bb=score(image,xa,ya); 
		before = (1.2*aa)+(bb); 
		swap (image,x,y,0,0,xa,ya);
		aa=score(image,x,y); 
		bb=score(image,xa,ya); 
		after = (aa)+(1.2*bb); 

		if (after>=before){
			swap (image,x,y,0,0,xa,ya);
			continue;
		} 

		imagel[x+(y*X_SIZE)]=30;
		imagel[xa+(ya*X_SIZE)]=30;



		if (f%1000000==0){
			disp(image,frame,1);
			frame++;
			printf ("%d \n",frame);
		}
		f++;


		//}}
	}

	scanf("%c",junk);


	close_x();
	exit(0);
}	

void disp (double *image2,int fram,int ab)
{
	int x,y;
	char *input;
	input=malloc(300);

	unsigned char *image;


        image=(unsigned char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer


       	for (y=0;y<Y_SIZE;y++)
       	{
               	int p=y*X_SIZE*3;
               	int XYP=X_SIZE*4*y;
               	for (x=0;x<X_SIZE;x++)
               	{
			int xpoint;
			int X_POINT;
			X_POINT=XYP+(4*x);
			xpoint=(x*3)+(p);

			x_buffer[X_POINT+2]=image2[xpoint];
			x_buffer[X_POINT+1]=image2[xpoint+1];
			x_buffer[X_POINT]=image2[xpoint+2];

			image[xpoint]=image2[xpoint];
			image[xpoint+1]=image2[xpoint+1];
			image[xpoint+2]=image2[xpoint+2];

                }
        }
	XPutImage(dis, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	sprintf(input,"./jpegs/rr%05d.jpg",fram);
	if (ab){jayit(image,X_SIZE, Y_SIZE, input);}
	free (input);
	free(image);
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
        unsigned long black,white;

        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //y_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //z_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        dis=XOpenDisplay((char *)0);
        screen=DefaultScreen(dis);
        black=BlackPixel(dis,screen),
        white=WhitePixel(dis,screen);
        win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        XSetStandardProperties(dis,win,"image","images",None,NULL,0,NULL);
        gc=XCreateGC(dis, win, 0,0);
        XSetBackground(dis,gc,black); XSetForeground(dis,gc,white);
        XClearWindow(dis, win);
        XMapRaised(dis, win);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(dis, 0);
        x_image=XCreateImage(dis, visual, DefaultDepth(dis,DefaultScreen(dis)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(dis, gc);
        XDestroyWindow(dis,win);
        XCloseDisplay(dis);
        exit(1);
};

void redraw() {
        XClearWindow(dis, win);
};

