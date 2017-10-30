#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
/* adaptive rejection metropolis sampling */

/* *********************************************************************** */

#include           <stdio.h>
#include           <math.h>
#include           <stdlib.h>

/* *********************************************************************** */

typedef struct point {    /* a point in the x,y plane */
double x,y;             /* x and y coordinates */
double ey;              /* exp(y-ymax+YCEIL) */
double cum;             /* integral up to x of rejection envelope */
int f;                  /* is y an evaluated point of log-density */
struct point *pl,*pr;   /* envelope points to left and right of x */
} POINT;

/* *********************************************************************** */

typedef struct envelope {  /* attributes of the entire rejection envelope */
int cpoint;              /* number of POINTs in current envelope */
int npoint;              /* max number of POINTs allowed in envelope */
int *neval;              /* number of function evaluations performed */
double ymax;             /* the maximum y-value in the current envelope */
POINT *p;                /* start of storage of envelope POINTs */
double *convex;          /* adjustment for convexity */
} ENVELOPE;

/* *********************************************************************** */

typedef struct funbag { /* everything for evaluating log density          */
void *mydata;      /* user-defined structure holding data for density */
double (*myfunc)(double x, void *mydata);
/* user-defined function evaluating log density at x */
} FUNBAG;

/* *********************************************************************** */

typedef struct metropolis { /* for metropolis step */
int on;            /* whether metropolis is to be used */
double xprev;      /* previous Markov chain iterate */
double yprev;      /* current log density at xprev */
} METROPOLIS;

/* *********************************************************************** */

//#define RAND_MAX 2147483647      /* For Sun4 : remove this for some systems */ 
#define XEPS  0.00001            /* critical relative x-value difference */
#define YEPS  0.1                /* critical y-value difference */
#define EYEPS 0.001              /* critical relative exp(y) difference */
#define YCEIL 50.                /* maximum y avoiding overflow in exp(y) */

/* *********************************************************************** */

/* declarations for functions defined in this file */

int arms_simple (int ninit, double *xl, double *xr, 
                 double (*myfunc)(double x, void *mydata), void *mydata,
                 int dometrop, double *xprev, double *xsamp);

int arms (double *xinit, int ninit, double *xl, double *xr, 
          double (*myfunc)(double x, void *mydata), void *mydata,
          double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
          int nsamp, double *qcent, double *xcent, int ncent,
          int *neval);

int initial (double *xinit, int ninit, double xl, double xr, int npoint,
             FUNBAG *lpdf, ENVELOPE *env, double *convex, int *neval,
             METROPOLIS *metrop);

void sample(ENVELOPE *env, POINT *p);

void invert(double prob, ENVELOPE *env, POINT *p);

int test(ENVELOPE *env, POINT *p, FUNBAG *lpdf, METROPOLIS *metrop);

int update(ENVELOPE *env, POINT *p, FUNBAG *lpdf, METROPOLIS *metrop);

void cumulate(ENVELOPE *env);

int meet (POINT *q, ENVELOPE *env, METROPOLIS *metrop);

double area(POINT *q);

double expshift(double y, double y0);

double logshift(double y, double y0);

double perfunc(FUNBAG *lpdf, ENVELOPE *env, double x);

void display(FILE *f, ENVELOPE *env);

double u_random();

/* *********************************************************************** */

int arms_simple (int ninit, double *xl, double *xr,
                 double (*myfunc)(double x, void *mydata), void *mydata,
                 int dometrop, double *xprev, double *xsamp)
  
  /* adaptive rejection metropolis sampling - simplified argument list */
  /* ninit        : number of starting values to be used */
  /* *xl          : left bound */
  /* *xr          : right bound */
  /* *myfunc      : function to evaluate log-density */
  /* *mydata      : data required by *myfunc */
  /* dometrop     : whether metropolis step is required */
  /* *xprev       : current value from markov chain */
  /* *xsamp       : to store sampled value */
  
{
  double xinit[ninit], convex=1.0, qcent, xcent;
  int err, i, npoint=100, nsamp=1, ncent=0, neval; 
  
  /* set up starting values */
  for(i=0; i<ninit; i++){
    xinit[i] = *xl + (i + 1.0) * (*xr - *xl)/(ninit + 1.0);
  }
  
  err = arms(xinit,ninit,xl,xr,myfunc,mydata,&convex,npoint,dometrop,xprev,xsamp,
             nsamp,&qcent,&xcent,ncent,&neval);
  
  return err;
}

/* *********************************************************************** */

int arms (double *xinit, int ninit, double *xl, double *xr, 
          double (*myfunc)(double x, void *mydata), void *mydata,
          double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
          int nsamp, double *qcent, double *xcent,
          int ncent, int *neval)
  
  /* to perform derivative-free adaptive rejection sampling with metropolis step */
  /* *xinit       : starting values for x in ascending order */
  /* ninit        : number of starting values supplied */
  /* *xl          : left bound */
  /* *xr          : right bound */
  /* *myfunc      : function to evaluate log-density */
  /* *mydata      : data required by *myfunc */
  /* *convex      : adjustment for convexity */
  /* npoint       : maximum number of envelope points */
  /* dometrop     : whether metropolis step is required */
  /* *xprev       : previous value from markov chain */
  /* *xsamp       : to store sampled values */
  /* nsamp        : number of sampled values to be obtained */
  /* *qcent       : percentages for envelope centiles */
  /* *xcent       : to store requested centiles */
  /* ncent        : number of centiles requested */
  /* *neval       : on exit, the number of function evaluations performed */
  
{
  
  ENVELOPE *env;      /* rejection envelope */
  POINT pwork;        /* a working point, not yet incorporated in envelope */
  int msamp=0;        /* the number of x-values currently sampled */
  FUNBAG lpdf;        /* to hold density function and its data */
  METROPOLIS *metrop; /* to hold bits for metropolis step */
  int i,err;
  
  /* check requested envelope centiles */
  for(i=0; i<ncent; i++){
    if((qcent[i] < 0.0) || (qcent[i] > 100.0)){
      /* percentage requesting centile is out of range */
      return 1005;
    }
  }
  
  /* incorporate density function and its data into FUNBAG lpdf */
  lpdf.mydata = mydata;
  lpdf.myfunc = myfunc;
  
  /* set up space required for envelope */
  env = (ENVELOPE *)malloc(sizeof(ENVELOPE));
  if(env == NULL){
    /* insufficient space */
    return 1006;
  }
  
  /* start setting up metropolis struct */
  metrop = (METROPOLIS *)malloc(sizeof(METROPOLIS));
  if(metrop == NULL){
    /* insufficient space */
    return 1006;
  }
  metrop->on = dometrop;
  
  /* set up initial envelope */
  err = initial(xinit,ninit,*xl,*xr,npoint,&lpdf,env,convex,
                neval,metrop);
  if(err)return err;
  
  /* finish setting up metropolis struct (can only do this after */
  /* setting up env) */
  if(metrop->on){
    if((*xprev < *xl) || (*xprev > *xr)){
      /* previous markov chain iterate out of range */
      return 1007;
    }
    metrop->xprev = *xprev;
    metrop->yprev = perfunc(&lpdf,env,*xprev);
  }
  
  /* now do adaptive rejection */
  do {
    /* sample a new point */
    sample (env,&pwork);
    
    /* perform rejection (and perhaps metropolis) tests */
    i = test(env,&pwork,&lpdf,metrop);
    if(i == 1){
      /* point accepted */
      xsamp[msamp++] = pwork.x;
    } else if (i != 0) {
      /* envelope error - violation without metropolis */
      return 2000;
    }  
  } while (msamp < nsamp);
  
  /* nsamp points now sampled */
  /* calculate requested envelope centiles */
  for (i=0; i<ncent; i++){
    invert(qcent[i]/100.0,env,&pwork);
    xcent[i] = pwork.x;
  }
  
  /* free space */
  free(env->p);
  free(env);
  free(metrop);
  
  return 0;
}

/* *********************************************************************** */

int initial (double *xinit, int ninit, double xl, double xr, int npoint,
             FUNBAG *lpdf, ENVELOPE *env, double *convex, int *neval,
             METROPOLIS *metrop)
  
  /* to set up initial envelope */
  /* xinit        : initial x-values */
  /* ninit        : number of initial x-values */
  /* xl,xr        : lower and upper x-bounds */
  /* npoint       : maximum number of POINTs allowed in envelope */
  /* *lpdf        : to evaluate log density */
  /* *env         : rejection envelope attributes */
  /* *convex      : adjustment for convexity */
  /* *neval       : current number of function evaluations */
  /* *metrop      : for metropolis step */
  
{
  int i,j,k,mpoint;
  POINT *q;
  
  if(ninit<3){
    /* too few initial points */
    return 1001;
  }
  
  mpoint = 2*ninit + 1;
  if(npoint < mpoint){
    /* too many initial points */
    return 1002;
  }
  
  if((xinit[0] <= xl) || (xinit[ninit-1] >= xr)){
    /* initial points do not satisfy bounds */
    return 1003;
  }
  
  for(i=1; i<ninit; i++){
    if(xinit[i] <= xinit[i-1]){
      /* data not ordered */
      return 1004;
    }
  }
  
  if(*convex < 0.0){
    /* negative convexity parameter */
    return 1008;
  }
  
  /* copy convexity address to env */
  env->convex = convex;
  
  /* copy address for current number of function evaluations */
  env->neval = neval;
  /* initialise current number of function evaluations */
  *(env->neval) = 0;
  
  /* set up space for envelope POINTs */
  env->npoint = npoint;
  env->p = (POINT *)malloc(npoint*sizeof(POINT));
  if(env->p == NULL){
    /* insufficient space */
    return 1006;
  }
  
  /* set up envelope POINTs */
  q = env->p;
  /* left bound */
  q->x = xl;
  q->f = 0;
  q->pl = NULL;
  q->pr = q+1;
  for(j=1, k=0; j<mpoint-1; j++){
    q++;
    if(j%2){
      /* point on log density */
      q->x = xinit[k++];
      q->y = perfunc(lpdf,env,q->x);
      q->f = 1;
    } else {
      /* intersection point */
      q->f = 0;
    }
    q->pl = q-1;
    q->pr = q+1;
  }
  /* right bound */
  q++;
  q->x = xr;
  q->f = 0;
  q->pl = q-1;
  q->pr = NULL;
  
  /* calculate intersection points */
  q = env->p;
  for (j=0; j<mpoint; j=j+2, q=q+2){
    if(meet(q,env,metrop)){
      /* envelope violation without metropolis */
      return 2000;
    }
  }
  
  /* exponentiate and integrate envelope */
  cumulate(env);
  
  /* note number of POINTs currently in envelope */
  env->cpoint = mpoint;
  
  return 0;
}

/* *********************************************************************** */

void sample(ENVELOPE *env, POINT *p)
  
  /* To sample from piecewise exponential envelope */
  /* *env    : envelope attributes */
  /* *p      : a working POINT to hold the sampled value */
  
{
  double prob;
  
  /* sample a uniform */
  prob = u_random();
  /* get x-value correponding to a cumulative probability prob */
  invert(prob,env,p);
  
  return;
}

/* *********************************************************************** */

void invert(double prob, ENVELOPE *env, POINT *p)
  
  /* to obtain a point corresponding to a qiven cumulative probability */
  /* prob    : cumulative probability under envelope */
  /* *env    : envelope attributes */
  /* *p      : a working POINT to hold the sampled value */
  
{
  double u,xl,xr,yl,yr,eyl,eyr,prop,z;
  POINT *q;
  
  /* find rightmost point in envelope */
  q = env->p;
  while(q->pr != NULL)q = q->pr;
  
  /* find exponential piece containing point implied by prob */
  u = prob * q->cum;
  while(q->pl->cum > u)q = q->pl;
  
  /* piece found: set left and right POINTs of p, etc. */
  p->pl = q->pl;
  p->pr = q;
  p->f = 0;
  p->cum = u;
  
  /* calculate proportion of way through integral within this piece */
  prop = (u - q->pl->cum) / (q->cum - q->pl->cum);
  
  /* get the required x-value */
  if (q->pl->x == q->x){
    /* interval is of zero length */
    p->x = q->x;
    p->y = q->y;
    p->ey = q->ey;
  } else {
    xl = q->pl->x;
    xr = q->x;
    yl = q->pl->y;
    yr = q->y;
    eyl = q->pl->ey;
    eyr = q->ey;
    if(fabs(yr - yl) < YEPS){
      /* linear approximation was used in integration in function cumulate */
      if(fabs(eyr - eyl) > EYEPS*fabs(eyr + eyl)){
        p->x = xl + ((xr - xl)/(eyr - eyl))
        * (-eyl + sqrt((1. - prop)*eyl*eyl + prop*eyr*eyr));
      } else {
        p->x = xl + (xr - xl)*prop;
      }
      p->ey = ((p->x - xl)/(xr - xl)) * (eyr - eyl) + eyl;
      p->y = logshift(p->ey, env->ymax);
    } else {
      /* piece was integrated exactly in function cumulate */
      p->x = xl + ((xr - xl)/(yr - yl))
      * (-yl + logshift(((1.-prop)*eyl + prop*eyr), env->ymax));
      p->y = ((p->x - xl)/(xr - xl)) * (yr - yl) + yl;
      p->ey = expshift(p->y, env->ymax);
    }
  }
  
  /* guard against imprecision yielding point outside interval */
  if ((p->x < xl) || (p->x > xr))exit(1);
  
  return;
}

/* *********************************************************************** */

int test(ENVELOPE *env, POINT *p, FUNBAG *lpdf, METROPOLIS *metrop)
  
  /* to perform rejection, squeezing, and metropolis tests */
  /* *env          : envelope attributes */
  /* *p            : point to be tested */
  /* *lpdf         : to evaluate log-density */
  /* *metrop       : data required for metropolis step */
  
{
  double u,y,ysqueez,ynew,yold,znew,zold,w;
  POINT *ql,*qr;
  
  /* for rejection test */
  u = u_random() * p->ey;
  y = logshift(u,env->ymax);
  
  if(!(metrop->on) && (p->pl->pl != NULL) && (p->pr->pr != NULL)){
    /* perform squeezing test */
    if(p->pl->f){
      ql = p->pl;
    } else {
      ql = p->pl->pl;
    }
    if(p->pr->f){
      qr = p->pr;
    } else {
      qr = p->pr->pr;
    }
    ysqueez = (qr->y * (p->x - ql->x) + ql->y * (qr->x - p->x))
      /(qr->x - ql->x);
    if(y <= ysqueez){
      /* accept point at squeezing step */
      return 1;
    }
  }
  
  /* evaluate log density at point to be tested */
  ynew = perfunc(lpdf,env,p->x);
  
  /* perform rejection test */
  if(!(metrop->on) || ((metrop->on) && (y >= ynew))){
    /* update envelope */
    p->y = ynew;
    p->ey = expshift(p->y,env->ymax);
    p->f = 1;
    if(update(env,p,lpdf,metrop)){
      /* envelope violation without metropolis */
      return -1;
    }
    /* perform rejection test */
    if(y >= ynew){
      /* reject point at rejection step */
      return 0;
    } else {
      /* accept point at rejection step */
      return 1;
    }
  }
  
  /* continue with metropolis step */
  yold = metrop->yprev;
  /* find envelope piece containing metrop->xprev */
  ql = env->p;
  while(ql->pl != NULL)ql = ql->pl;
  while(ql->pr->x < metrop->xprev)ql = ql->pr;
  qr = ql->pr;
  /* calculate height of envelope at metrop->xprev */
  w = (metrop->xprev - ql->x)/(qr->x - ql->x);
  zold = ql->y + w*(qr->y - ql->y);
  znew = p->y;
  if(yold < zold)zold = yold;
  if(ynew < znew)znew = ynew;
  w = ynew-znew-yold+zold;
  if(w > 0.0)w = 0.0;
  
  if(w > -YCEIL){
    w = exp(w);
  } else {
    w = 0.0;
  }
  u = u_random();
  if(u > w){
    /* metropolis says dont move, so replace current point with previous */
    /* markov chain iterate */
    p->x = metrop->xprev;
    p->y = metrop->yprev;
    p->ey = expshift(p->y,env->ymax);
    p->f = 1;
    p->pl = ql;
    p->pr = qr;
  } else {
    /* trial point accepted by metropolis, so update previous markov */
    /* chain iterate */
    metrop->xprev = p->x;
    metrop->yprev = ynew;
  }
  return 1;
}

/* *********************************************************************** */

int update(ENVELOPE *env, POINT *p, FUNBAG *lpdf, METROPOLIS *metrop)
  
  /* to update envelope to incorporate new point on log density*/
  /* *env          : envelope attributes */
  /* *p            : point to be incorporated */
  /* *lpdf         : to evaluate log-density */
  /* *metrop       : for metropolis step */
  
{
  POINT *m,*ql,*qr,*q;
  
  if(!(p->f) || (env->cpoint > env->npoint - 2)){
    /* y-value has not been evaluated or no room for further points */
    /* ignore this point */
    return 0;
  }
  
  /* copy working POINT p to a new POINT q */
  q = env->p + env->cpoint++;
  q->x = p->x;
  q->y = p->y;
  q->f = 1;
  
  /* allocate an unused POINT for a new intersection */
  m = env->p + env->cpoint++;
  m->f = 0;
  if((p->pl->f) && !(p->pr->f)){
    /* left end of piece is on log density; right end is not */
    /* set up new intersection in interval between p->pl and p */
    m->pl = p->pl;
    m->pr = q;
    q->pl = m;
    q->pr = p->pr;
    m->pl->pr = m;
    q->pr->pl = q;
  } else if (!(p->pl->f) && (p->pr->f)){
    /* left end of interval is not on log density; right end is */
    /* set up new intersection in interval between p and p->pr */
    m->pr = p->pr;
    m->pl = q;
    q->pr = m;
    q->pl = p->pl;
    m->pr->pl = m;
    q->pl->pr = q;
  } else {
    /* this should be impossible */
    exit(10);
  }
  
  /* now adjust position of q within interval if too close to an endpoint */
  if(q->pl->pl != NULL){
    ql = q->pl->pl;
  } else {
    ql = q->pl;
  }
  if(q->pr->pr != NULL){
    qr = q->pr->pr;
  } else {
    qr = q->pr;
  }
  if (q->x < (1. - XEPS) * ql->x + XEPS * qr->x){
    /* q too close to left end of interval */
    q->x = (1. - XEPS) * ql->x + XEPS * qr->x;
    q->y = perfunc(lpdf,env,q->x);
  } else if (q->x > XEPS * ql->x + (1. - XEPS) * qr->x){
    /* q too close to right end of interval */
    q->x = XEPS * ql->x + (1. - XEPS) * qr->x;
    q->y = perfunc(lpdf,env,q->x);
  }
  
  /* revise intersection points */
  if(meet(q->pl,env,metrop)){
    /* envelope violation without metropolis */
    return 1;
  }
  if(meet(q->pr,env,metrop)){
    /* envelope violation without metropolis */
    return 1;
  }
  if(q->pl->pl != NULL){
    if(meet(q->pl->pl->pl,env,metrop)){
      /* envelope violation without metropolis */
      return 1;
    }
  }
  if(q->pr->pr != NULL){
    if(meet(q->pr->pr->pr,env,metrop)){
      /* envelope violation without metropolis */
      return 1;
    }
  }
  
  /* exponentiate and integrate new envelope */
  cumulate(env);
  
  return 0;
}

/* *********************************************************************** */

void cumulate(ENVELOPE *env)
  
  /* to exponentiate and integrate envelope */
  /* *env     : envelope attributes */
  
{
  POINT *q,*qlmost;
  
  qlmost = env->p;
  /* find left end of envelope */
  while(qlmost->pl != NULL)qlmost = qlmost->pl;
  
  /* find maximum y-value: search envelope */
  env->ymax = qlmost->y;
  for(q = qlmost->pr; q != NULL; q = q->pr){
    if(q->y > env->ymax)env->ymax = q->y;
  }
  
  /* exponentiate envelope */
  for(q = qlmost; q != NULL; q = q->pr){
    q->ey = expshift(q->y,env->ymax);
  }
  
  /* integrate exponentiated envelope */
  qlmost->cum = 0.;
  for(q = qlmost->pr; q != NULL; q = q->pr){
    q->cum = q->pl->cum + area(q);
  }
  
  return;
}

/* *********************************************************************** */

int meet (POINT *q, ENVELOPE *env, METROPOLIS *metrop)
  /* To find where two chords intersect */
  /* q         : to store point of intersection */
  /* *env      : envelope attributes */
  /* *metrop   : for metropolis step */
  
{
  double gl,gr,grl,dl,dr;
  int il,ir,irl;
  
  if(q->f){
    /* this is not an intersection point */
    exit(30);
  }
  
  /* calculate coordinates of point of intersection */
  if ((q->pl != NULL) && (q->pl->pl->pl != NULL)){
    /* chord gradient can be calculated at left end of interval */
    gl = (q->pl->y - q->pl->pl->pl->y)/(q->pl->x - q->pl->pl->pl->x);
    il = 1;
  } else {
    /* no chord gradient on left */
    il = 0;
  }
  if ((q->pr != NULL) && (q->pr->pr->pr != NULL)){
    /* chord gradient can be calculated at right end of interval */
    gr = (q->pr->y - q->pr->pr->pr->y)/(q->pr->x - q->pr->pr->pr->x);
    ir = 1;
  } else {
    /* no chord gradient on right */
    ir = 0;
  }
  if ((q->pl != NULL) && (q->pr != NULL)){
    /* chord gradient can be calculated across interval */
    grl = (q->pr->y - q->pl->y)/(q->pr->x - q->pl->x);
    irl = 1;
  } else {
    irl = 0;
  }
  
  if(irl && il && (gl<grl)){
    /* convexity on left exceeds current threshold */
    if(!(metrop->on)){
      /* envelope violation without metropolis */
      return 1;
    }
    /* adjust left gradient */
    gl = gl + (1.0 + *(env->convex)) * (grl - gl);
  }
  
  if(irl && ir && (gr>grl)){
    /* convexity on right exceeds current threshold */
    if(!(metrop->on)){
      /* envelope violation without metropolis */
      return 1;
    }
    /* adjust right gradient */
    gr = gr + (1.0 + *(env->convex)) * (grl - gr);
  }
  
  if(il && irl){
    dr = (gl - grl) * (q->pr->x - q->pl->x);
    if(dr < YEPS){
      /* adjust dr to avoid numerical problems */
      dr = YEPS;
    }
  }
  
  if(ir && irl){
    dl = (grl - gr) * (q->pr->x - q->pl->x);
    if(dl < YEPS){
      /* adjust dl to avoid numerical problems */
      dl = YEPS;
    }
  }
  
  if(il && ir && irl){
    /* gradients on both sides */
    q->x = (dl * q->pr->x + dr * q->pl->x)/(dl + dr);
    q->y = (dl * q->pr->y + dr * q->pl->y + dl * dr)/(dl + dr);
  } else if (il && irl){
    /* gradient only on left side, but not right hand bound */
    q->x = q->pr->x;
    q->y = q->pr->y + dr;
  } else if (ir && irl){
    /* gradient only on right side, but not left hand bound */
    q->x = q->pl->x;
    q->y = q->pl->y + dl;
  } else if (il){
    /* right hand bound */
    q->y = q->pl->y + gl * (q->x - q->pl->x);
  } else if (ir){
    /* left hand bound */
    q->y = q->pr->y - gr * (q->pr->x - q->x);
  } else {
    /* gradient on neither side - should be impossible */
    exit(31);
  }
  if(((q->pl != NULL) && (q->x < q->pl->x)) ||
     ((q->pr != NULL) && (q->x > q->pr->x))){
    /* intersection point outside interval (through imprecision) */
    exit(32);
  }
  /* successful exit : intersection has been calculated */
  return 0;
}

/* *********************************************************************** */

double area(POINT *q)
  
  /* To integrate piece of exponentiated envelope to left of POINT q */
  
{
  double a;
  
  if(q->pl == NULL){
    /* this is leftmost point in envelope */
    exit(1);
  } else if(q->pl->x == q->x){
    /* interval is zero length */
    a = 0.;
  } else if (fabs(q->y - q->pl->y) < YEPS){
    /* integrate straight line piece */
    a = 0.5*(q->ey + q->pl->ey)*(q->x - q->pl->x);
  } else {
    /* integrate exponential piece */
    a = ((q->ey - q->pl->ey)/(q->y - q->pl->y))*(q->x - q->pl->x);
  }
  return a;
}

/* *********************************************************************** */

double expshift(double y, double y0)
  
  /* to exponentiate shifted y without underflow */
{
  if(y - y0 > -2.0 * YCEIL){
    return exp(y - y0 + YCEIL);
  } else {
    return 0.0;
  }
}

/* *********************************************************************** */

double logshift(double y, double y0)
  
  /* inverse of function expshift */
{
  return (log(y) + y0 - YCEIL);
}

/* *********************************************************************** */

double perfunc(FUNBAG *lpdf, ENVELOPE *env, double x)
  
  /* to evaluate log density and increment count of evaluations */
  
  /* *lpdf   : structure containing pointers to log-density function and data */
  /* *env    : envelope attributes */
  /* x       : point at which to evaluate log density */
  
{
  double y;
  
  /* evaluate density function */
  y = (lpdf->myfunc)(x,lpdf->mydata);
  
  /* increment count of function evaluations */
  (*(env->neval))++;
  
  return y;
}

/* *********************************************************************** */

// void display(FILE *f, ENVELOPE *env)
//   
//   /* to display envelope - for debugging only */
// {
//   POINT *q;
//   
//   /* print envelope attributes */
//   fprintf(f,"========================================================\n");
//   fprintf(f,"envelope attributes:\n");
//   fprintf(f,"points in use = %d, points available = %d\n",
//           env->cpoint,env->npoint);
//   fprintf(f,"function evaluations = %d\n",*(env->neval));
//   fprintf(f,"ymax = %f, p = %x\n",env->ymax,env->p);
//   fprintf(f,"convexity adjustment = %f\n",*(env->convex));
//   fprintf(f,"--------------------------------------------------------\n");
//   
//   /* find leftmost POINT */
//   q = env->p;
//   while(q->pl != NULL)q = q->pl;
//   
//   /* now print each POINT from left to right */
//   for(q = env->p; q != NULL; q = q->pr){
//     fprintf(f,"point at %x, left at %x, right at %x\n",q,q->pl,q->pr);
//     fprintf(f,"x = %f, y = %f, ey = %f, cum = %f, f = %d\n",
//             q->x,q->y,q->ey,q->cum,q->f);
//   }
//   fprintf(f,"========================================================\n");
//   
//   return;
// }

/* *********************************************************************** */

double u_random()
  
  /* to return a standard uniform random number */
{
  return ((double)rand() + 0.5)/((double)RAND_MAX + 1.0);
}

/* *********************************************************************** */
struct norm_parm {
  double mean,sd;
};
struct params {
  double  alpha, beta, s_t, u_j, sigma2,  y_tj, n_tj;
};
/* ********************************************************************* */
/*Sampling V*/
/* ********************************************************************* */
struct v_para{
  double y_ij, n_ij, z_ij, delta0;
};
double vij_ld(double v_ij, void *v_para_input){
  struct v_para *d = static_cast<v_para*>(v_para_input);
  double ld = v_ij * d->y_ij - d->n_ij * log(1+exp(v_ij)) - (0.5/d->delta0) * pow(v_ij - d->z_ij, 2);
  return ld;
};
// [[Rcpp::export]]
arma::mat sample_v(arma::mat v_mat, arma::mat y_mat, arma::mat n_mat, arma::mat z_mat, double delta0,
                   double xl, double xr, int I, int J){
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
  int neval;
  double xinit[ninit];
  for(int i=0;i<ninit;i++){
    xinit[i] = xl + (i + 1.0) * (xr - xl)/(ninit + 1.0);
  }
  double xsamp[100], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  int dometrop = 1;
  
  arma::mat v_mat_new = arma::zeros<arma::mat>(I, J);
  for(int i=0; i<I; i++){
    for(int j=0; j<J; j++){
      struct v_para v_para_input;
      v_para_input.delta0 = delta0;
      v_para_input.n_ij = n_mat(i,j);
      v_para_input.y_ij = y_mat(i,j);
      v_para_input.z_ij = z_mat(i,j);
      double xprev = v_mat(i,j);
      err = arms(xinit,ninit,&xl,&xr,vij_ld,&v_para_input,&convex,
                 npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
      v_mat_new(i,j) = xsamp[0];
    }
  }
  return v_mat_new;
};
/* ********************************************************************* */
/*Sampling z */
/* ********************************************************************* */
// [[Rcpp::export]]
arma::vec sample_z(arma::vec v_vec, double delta0, double rho, arma::mat Sigma_inv, arma::mat C, int I, int J){
  arma::mat B = arma::eye<arma::mat>(I,I) - rho *C;
  arma::mat H = arma::kron(Sigma_inv, B);
  arma::mat I_delta0_H = arma::eye<arma::mat>(I*J,I*J) + delta0 * H;
  arma::mat I_delta0_H_inv = arma::inv_sympd(I_delta0_H);
  arma::mat Cov = delta0 * I_delta0_H_inv;
  arma::vec mu_star = I_delta0_H_inv*v_vec;
  //mvnorm
  arma::mat R = arma::chol(Cov);
  arma::vec z_new = R.t()*arma::randn(I*J) + mu_star;
  return z_new;
}
/* ********************************************************************* */
/*Sampling r */
/* ********************************************************************* */
// [[Rcpp::export]]
arma::vec sample_r(arma::vec v_vec, double delta0, double rho, arma::mat Sigma_inv, arma::mat C, int I, int J,
                   int p, arma::mat X, double delta_beta){
  int IJ = I*J;
  arma::mat B = arma::eye<arma::mat>(I,I) - rho *C;
  arma::mat H = arma::kron(Sigma_inv, B);
  arma::mat I_delta0_H = arma::eye<arma::mat>(I*J,I*J) + delta0 * H;
  arma::mat Omega_r = arma::zeros<arma::mat>(IJ+p, IJ+p);
  Omega_r.submat(0,0,p-1,p-1) = 1.0/delta_beta * arma::eye<arma::mat>(p,p);
  Omega_r.submat(p,p, IJ+p-1, IJ+p-1) = H;
  arma::mat Omega_p_part = delta0 * Omega_r + X.t()*X;
  arma::mat Sigma_p_part = arma::inv_sympd(Omega_p_part);
  arma::mat Sigma_p = delta0 * Sigma_p_part;
  arma::vec mu_p = Sigma_p_part * X.t() * v_vec;
  //mvnorm
  arma::mat R = arma::chol(Sigma_p);
  arma::vec r_new = R.t()*arma::randn(IJ+p) + mu_p;
  return r_new;
}
/* ********************************************************************* */
/*Sampling r_sp */
/* ********************************************************************* */
// [[Rcpp::export]]
arma::vec sample_r_sp(arma::vec v_vec, double delta0, double rho, arma::mat C, int I, int J,
                      int p, arma::mat X, double delta_beta, double eta1){
  arma::mat beta_part = delta_beta*arma::eye<arma::mat>(p,p);
  double delta1 = delta0/eta1;
  arma::mat IrhoC = arma::eye<arma::mat>(I,I) - rho * C;
  arma::mat inv_IrhoC = arma::inv_sympd(IrhoC);
  arma::mat w_part = delta1 * inv_IrhoC;
  arma::mat Sigma = arma::zeros<arma::mat>(I+p, I+p);
  Sigma.submat(0,0,p-1,p-1) = beta_part;
  Sigma.submat(p,p, I+p-1, I+p-1) = w_part;
  arma::mat Sigma_inv = arma::zeros<arma::mat>(I+p, I+p);
  Sigma_inv.submat(0,0,p-1,p-1) = (1.0/delta_beta)*arma::eye<arma::mat>(p,p);
  Sigma_inv.submat(p,p, I+p-1, I+p-1) = (1.0/delta1)*IrhoC;
  
  arma::mat Sigma_star = arma::inv_sympd(delta0*Sigma_inv + X.t()*X);
  arma::mat Cov = delta0 * Sigma_star;
  arma::vec mu_star = Sigma_star*X.t()*v_vec;
  //mvnorm
  arma::mat R = arma::chol(Cov);
  arma::vec r_new = R.t()*arma::randn(I+p) + mu_star;
  return r_new;
}

/* ********************************************************************* */
/*Sampling rho */
/* ********************************************************************* */
struct rho_para{
  arma::vec w;//the spatial vector;
  arma::mat C;
  arma::vec nu;
  double eta1, delta0;
  int K;
};

double rho_ld(double rho, void *rho_para_input){
  struct rho_para *d = static_cast<rho_para*>(rho_para_input);
  arma::vec det_vec = arma::ones<arma::vec>(d->K) - rho*d->nu; 
  double det_part = 0.5*arma::sum(arma::log(det_vec));
  double cross_prod = arma::as_scalar((d->w).t() * d->C * (d->w));
  double ld = 0.5 * rho * d->eta1/d->delta0 * cross_prod + det_part;
  return ld;
};
// [[Rcpp::export]]
double sample_rho(double xprev, arma::vec w_vec, arma::mat C, arma::vec nu_vec, 
                  double eta1, double delta0, int K, double rho_min, double rho_max){
  double xl=rho_min;
  double xr=rho_max;
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
  int neval;
  double xinit[ninit];
  for(int i=0;i<ninit;i++){
    xinit[i] = xl + (i + 1.0) * (xr - xl)/(ninit + 1.0);
  }
  double xsamp[100], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  int dometrop = 1;
  
  struct rho_para rho_para_input;
  rho_para_input.C = C;
  rho_para_input.delta0 = delta0;
  rho_para_input.eta1 = eta1;
  rho_para_input.K = K;
  rho_para_input.nu = nu_vec;
  rho_para_input.w = w_vec;
  
  err = arms(xinit,ninit,&xl,&xr,rho_ld,&rho_para_input,&convex,
             npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
  if(err>0){
    exit(1);
  }
  double rho_new = xsamp[0];
  return rho_new;

}
/* ********************************************************************* */
/*Sampling eta1 */
/* ********************************************************************* */
struct eta1_para{
  arma::vec w;
  arma::mat C;
  double delta0, rho;
  int K;
};
//sample phi=log(eta1)
double phi_ld(double phi, void *eta1_para_input){
  struct eta1_para *d = static_cast<eta1_para*>(eta1_para_input);
  arma::mat IrhoC = arma::eye<arma::mat>(d->K,d->K) - d->rho * d->C;
  double WICW = arma::as_scalar((d->w).t()*IrhoC*(d->w));
  double ephi = exp(phi);
  double ld = phi*(0.5*d->K + 1) - 2*log(ephi+1) - 0.5*ephi/d->delta0 * WICW;
  return ld;
};
double eta1_ld(double eta1, void *eta1_para_input){
  struct eta1_para *d = static_cast<eta1_para*>(eta1_para_input);
  arma::mat IrhoC = arma::eye<arma::mat>(d->K,d->K) - d->rho * d->C;
  double WICW = arma::as_scalar((d->w).t()*IrhoC*(d->w));
  double ld = 0.5*d->K*log(eta1) - 2.0*log(eta1+1) - 0.5*eta1/d->delta0 * WICW;
  return ld;
}
// [[Rcpp::export]]
double sample_eta1(double xprev, arma::vec w_vec, arma::mat C, 
                   double rho, double delta0, int K, double xl, double xr){
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
  int neval;
  double xinit[ninit];
  for(int i=0;i<ninit;i++){
    xinit[i] = xl + (i + 1.0) * (xr - xl)/(ninit + 1.0);
  }
  double xsamp[100], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  int dometrop = 1;
  
  struct eta1_para eta1_para_input;
  eta1_para_input.w = w_vec;
  eta1_para_input.C = C;
  eta1_para_input.delta0 = delta0;
  eta1_para_input.rho = rho;
  eta1_para_input.K = K;
  err = arms(xinit,ninit,&xl,&xr,eta1_ld,&eta1_para_input,&convex,
             npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
  if(err>0){
    exit(1);
  }
  double eta1_new = xsamp[0];
  return eta1_new;
};
// [[Rcpp::export]]
double sample_eta1_by_log(double xprev, arma::vec w_vec, arma::mat C, 
                          double rho, double delta0, int K, double xl, double xr){
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
  int neval;
  double xinit[ninit];
  for(int i=0;i<ninit;i++){
    xinit[i] = xl + (i + 1.0) * (xr - xl)/(ninit + 1.0);
  }
  double xsamp[100], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  int dometrop = 1;
  double log_xprev = log(xprev);
  
  struct eta1_para eta1_para_input;
  eta1_para_input.w = w_vec;
  eta1_para_input.C = C;
  eta1_para_input.delta0 = delta0;
  eta1_para_input.rho = rho;
  eta1_para_input.K = K;
  err = arms(xinit,ninit,&xl,&xr,phi_ld,&eta1_para_input,&convex,
             npoint,dometrop,&log_xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
  if(err>0){
    exit(1);
  }
  double phi_new = xsamp[0];
  double eta1_new = exp(phi_new);
  return eta1_new;
};






