/*

  progressbar: A simple progress bar with elasped and estimated time
  to completion.

*/

#include<config.h>

#ifdef HAVE_TIME_H
 #include<time.h>
#endif

#ifdef HAVE_MATH_H
  #include<math.h>
#endif

#ifdef HAVE_STDIO_H
  #include<stdio.h>
#endif

int glob_nevents;
int glob_eventctr = 0;
int glob_curpos = 0;
int glob_maxpos = 25;
int glob_digits_x;

time_t glob_startbar_t;

void init_progressbar(int nevents)

/* We will move the bar each time an event is registered. The task
   should be completed after nevents. */

{
  int i;

  /* Draw headers for %complete, elapsed, estimated */

  printf("Progress                     %%comp  elapsed   remaining\n");

  /* Draw bar "blank" */

  printf("|");
  for(i=1;i<=glob_maxpos;i++) {
    printf("-");
  }
  printf("|");

  /* Now finish initializing globals */

  glob_curpos = 0;
  glob_startbar_t = time(NULL);
  glob_nevents = nevents;
  glob_eventctr = 0;

  /* Display */

  fflush(stdout);

}
  
void convert_to_hms(double dsec,int *hrs,int *min, int *sec)

{
  *hrs = (int)(floor(dsec/3600.0));
  dsec -= (*hrs)*3600;

  *min = (int)(floor(dsec/60.0));
  dsec -= (*min)*60;

  *sec  = (int)dsec;
}

void update_progressbar()

/* Updates the progress bar using globals. */

{
  int pcomp;
  time_t ctime;
  int newpos;
  int ehrs,emin,esec;
  int rhrs,rmin,rsec;
  double elapsed,remain;
  int i;

  glob_eventctr++;

  /* Update the bar itself. */

  printf("\r"); /* This returns to beginning of current line. Slick, huh? */

  newpos = (glob_maxpos * glob_eventctr)/glob_nevents;  
  
  printf("|");
  for (i=1;i<= newpos;i++) {
    printf("*");
  }
  for(;i<=glob_maxpos;i++) {
    printf("-");
  }
  printf("|");
  
  /* Update the numerical display */

  pcomp = (int)(100.0*(double)(glob_eventctr)/((double) glob_nevents));
  
  ctime = time(NULL);
  elapsed = difftime(ctime,glob_startbar_t);
  convert_to_hms(elapsed,&ehrs,&emin,&esec);
  remain = elapsed*((double)(glob_nevents)/(double)(glob_eventctr)) - elapsed;
  convert_to_hms(remain,&rhrs,&rmin,&rsec);

  printf("  %3d%%  %2d:%02d:%02d   %2d:%02d:%02d",
	   pcomp,ehrs,emin,esec,rhrs,rmin,rsec);

  /* And make sure it displays */

  fflush(stdout);

}
  
void free_progressbar()

{
  printf("\n"); /* Finally leave the current line */
}
  
  
  
  
