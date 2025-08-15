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

	sr=((avgr-image[ul])*(avgr-image[ul]))+((avgr-image[ur])*(avgr-image[ur]))+((avgr-image[um])*(avgr-image[um]))+
	   ((avgr-image[ml])*(avgr-image[ml]))+((avgr-image[mr])*(avgr-image[mr]))+((avgr-image[mm])*(avgr-image[mm]))+
	   ((avgr-image[dl])*(avgr-image[dl]))+((avgr-image[dr])*(avgr-image[dr]))+((avgr-image[dm])*(avgr-image[dm]));

	sg=((avgg-image[ul+1])*(avgg-image[ul+1]))+((avgg-image[ur+1])*(avgg-image[ur+1]))+((avgg-image[um+1])*(avgg-image[um+1]))+
	   ((avgg-image[ml+1])*(avgg-image[ml+1]))+((avgg-image[mr+1])*(avgg-image[mr+1]))+((avgg-image[mm+1])*(avgg-image[mm+1]))+
	   ((avgg-image[dl+1])*(avgg-image[dl+1]))+((avgg-image[dr+1])*(avgg-image[dr+1]))+((avgg-image[dm+1])*(avgg-image[dm+1]));


	sb=((avgb-image[ul+2])*(avgb-image[ul+2]))+((avgb-image[ur+2])*(avgb-image[ur+2]))+((avgb-image[um+2])*(avgb-image[um+2]))+
	   ((avgb-image[ml+2])*(avgb-image[ml+2]))+((avgb-image[mr+2])*(avgb-image[mr+2]))+((avgb-image[mm+2])*(avgb-image[mm+2]))+
	   ((avgb-image[dl+2])*(avgb-image[dl+2]))+((avgb-image[dr+2])*(avgb-image[dr+2]))+((avgb-image[dm+2])*(avgb-image[dm+2]));

	sd=sr+sg+sb;

	return sd;

}



void swap (double *image,int x, int y, int dxa, int dya, int dxb, int dyb)
{
	double sr,sg,sb;
	double dr,dg,db;

	sr=image[((x+dxa)*3)+(3*(y+dya)*X_SIZE)];
	sg=image[((x+dxa)*3)+(3*(y+dya)*X_SIZE)+1];
	sb=image[((x+dxa)*3)+(3*(y+dya)*X_SIZE)+2];

	dr=image[((x+dxb)*3)+(3*(y+dyb)*X_SIZE)];
	dg=image[((x+dxb)*3)+(3*(y+dyb)*X_SIZE)+1];
	db=image[((x+dxb)*3)+(3*(y+dyb)*X_SIZE)+2];

	image[((x+dxa)*3)+(3*(y+dya)*X_SIZE)]=dr;
	image[((x+dxa)*3)+(3*(y+dya)*X_SIZE)+1]=dg;
	image[((x+dxa)*3)+(3*(y+dya)*X_SIZE)+2]=db;

	image[((x+dxb)*3)+(3*(y+dyb)*X_SIZE)]=sr;
	image[((x+dxb)*3)+(3*(y+dyb)*X_SIZE)+1]=sg;
	image[((x+dxb)*3)+(3*(y+dyb)*X_SIZE)+2]=sb;

}

double reset (double *notes,int a, double top)
{

	double bottom;
	int i;

	//top=10000;
	bottom=8;

	if (a==-1){

	for (i=0;i<6;i++){ 
		notes[i]=bottom+((double)(rand()%(10000*(int)(top-bottom)))/10000);
	}} else{
		return bottom+((double)(rand()%(10000*(int)(top-bottom)))/10000);
	}
}

int main(int argc,char *argv[])
{
	double *image;
	double *wav;
	short *wavv;
	int i,j,loop,frame,bottom;
	double notes[6];
	double notel[6];
	int notef[3000];

	int rate;
	int dur;

	rate=48000;
	dur=60;


        image=(double *)malloc(sizeof (double)*X_SIZE*Y_SIZE*3); // disp buffer
        wav=(double *)malloc(sizeof (double)*2*dur*60*rate); // disp buffer
        wavv=(short *)malloc(sizeof (short)*2*dur*60*rate); // disp buffer


	char junk[30];


	init_x();



	int along,count,go,barren;
	double score,safe,best_score;
	long  wavp;
	along=0;
	best_score=0;
	count=0;
	barren=0;
	wavp=0;

	double angl[6],angr[6];

	for (i=0;i<6;i++){angr[i];angl[i]=0;notel[i]=0;}
	for (i=0;i<rate*dur*60*2;i++){wav[i]=0;}
	reset (notes,-1,200);


	int goodones,p;

	goodones=0;

	double top;

	top=200;
	

	while (goodones<150)
	{
		safe=notes[along];
		double sug;
		top+=0.003;

		sug=reset(notes,along,top);
		for (j=0;j<6;j++)
		{
			if (notes[j]==sug) { break;}
		}
		notes[along]=sug;
		score=0;
		for (i=0;i<6;i++)
		{
			if(i==along){continue;}
			double fit;
			if (notes[i]>sug){ 
				fit=(notes[i]/sug);
			}else{
				fit=(sug/notes[i]);
			}
			fit=(0.5-(fit-((int)fit)));
			//if (fit>0){ score+=fit;} else{score-=fit;}
			score+=fit*fit;
		}


		if (score>best_score){ 
			// we have a solution
			best_score=score;
			printf ("NS %f %d ",best_score,count);
			for (i=0;i<6;i++){
				printf (",%f",notes[i]);}
			printf ("\n");
			barren=0;


                        double amp;
                        for (i=0;i<6;i++){angl[i]=0;angr[i];}

                        amp=1;

			double prog;

			prog=(double)wavp/(2*60*60*rate);
                        p=rate/2+(rand()%(int)(rate*6*(1+prog)));

                        for (j=0;j<p;j++)
                        {
                                double left,right,last,next,il;
                                left=0;right=0;
				int m;
                                m=1+(2*(double)wavp/((double)rate*60*60));
				

                                next=(double)j/(double)(p);
                                last=1-next;
                                il=1-(prog);
                                //il=(prog);

                                for (i=0;i<6;i++)
                                {
                                        double l,r,lp,rp;
                                        double freq,lpan,rpan;

                                        lp=(double)(i)/5;
                                        rp=1-lp;

					lpan=(lp*next)+(rp*last);
					rpan=1-lpan;

                                        freq=(notes[i]);
                                        angl[i]+=2*M_PI*freq/(double)rate;
                                        angr[i]+=2*M_PI*freq*(int)(1+(prog*9))/(4*(double)rate);
                                        if (angl[i]>(2*M_PI)){ angl[i]-=2*M_PI;}
                                        if (angr[i]>(2*M_PI)){ angr[i]-=2*M_PI;}
                                        l=sin(angl[i]); r=sin(angr[i]);
                                        if (l>il){l=1;} if (l<-il){l=-1;}
                                        if (r>il){r=1;} if (r<-il){r=-1;}
                                        left+=(lpan*2276*l*amp)+((rpan)*2276*r*amp);
                                        right+=(rpan*2276*l*amp)+(lpan*2276*r*amp) ;
                                }
                                amp*=(0.999982 +(0.000009*(prog)));

                                wav[wavp]+=left;
                                wav[wavp+1]+=right;

				if (wavp<(60*60*2*rate)-(rate*120)){

                                wav[wavp+(int)(rate*9*prog)+2500]+=(7-(5*prog))*wav[wavp+1]/10;
                                wav[wavp+(int)(rate*11*prog)+2501]+=(7-(5*prog))*wav[wavp]/10; }
                                wavp+=2;
				if (wavp>=(60*60*2*rate)){goodones=1000;break;}
			}

		}else{ notes[along]=safe;}

		along++;
		if (along>5){along=0;}
		//scanf("%c",junk);
		barren++;
		count++;
		if (barren>10000){ 
			/*
			double amp;
			for (i=0;i<6;i++){ang[i]=0;}
			p=1+(rand()%(21));

			amp=1;
		        for (j=0;j<rate*p;j++)
                        {
                                double left,right,last,next,il;
                                left=0;right=0;

                                next=(double)j/(double)(rate*p);
                                last=1-next;
				il=1-((double)wavp/(rate*60*60));

                                for (i=0;i<6;i++)
                                {
                                        double l,r;
                                        double freq,pan;

					pan=(double)(i+1)/6;

					pan=(next*pan)+(last*(1-pan));

                                        freq=(notes[i]);
                                        ang[i]+=2*M_PI*freq/(double)rate;
                                        if (ang[i]>(2*M_PI)){ ang[i]-=2*M_PI;}
                                        l=sin(ang[i]); r=sin(ang[i]*3.01);
                                        if (l>il){l=1;} if (l<-il){l=-1;}
                                        if (r>il){r=1;} if (r<-il){r=-1;}
                                        left+=pan*5276*l*amp;
                                        right+=(1-pan)*5276*r*amp;
                                }
				amp*=(0.99999 +(0.0000099*((double)p/20)));
                                wav[wavp]+=left;
                                wav[wavp+1]+=right;

				wav[wavp+(rate*10)]+=9*wav[wavp+1]/10;
				wav[wavp+(rate/10)+1]+=9*wav[wavp]/10;
                                wavp+=2;

                        } */
                        for (i=0;i<6;i++){ notel[i]=notes[i];}
			
			barren=0;reset(notes,-1,top);count=0;score=0;along=0;goodones++;best_score=0;}
	}

	for (i=0;i<wavp;i++){wavv[i]=wav[i];}

	save_wav(wavv,"phase.wav", rate, 2, wavp );



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
	sprintf(input,"./jpegs/sf%05d.jpg",fram);
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

