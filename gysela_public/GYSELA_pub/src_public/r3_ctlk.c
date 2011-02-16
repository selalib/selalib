#include <fcntl.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/procfs.h>

#if defined (SGI_MIPS)
#include <sys/syssgi.h>
#include <sys/hwperftypes.h>
#include <sys/hwperfmacros.h>
#endif

#if defined (USE_R3_INFO_OMP)
#include <omp.h>
#endif
#if defined (USE_R3_INFO_MPI) || defined (USE_R3_INFO_SMA)
#include <mpi.h>
#endif

int  r3_fd = 0, r3_perf = 0;
char pfile[256];

double r3_secondr_(void)
{
   struct timespec t;

#if defined (USE_R3_INFO_OMP)
      return omp_get_wtime();
#elif defined (USE_R3_INFO_MPI) || defined (USE_R3_INFO_SMA)
      return MPI_Wtime();
#else

#if defined (SGI_MIPS)
      clock_gettime(CLOCK_SGI_CYCLE, &t);
#else
      clock_gettime(CLOCK_REALTIME, &t);
#endif

      return ((double)t.tv_sec+(double)t.tv_nsec*1.e-9);
#endif

}


#if defined (SGI_MIPS)
double r3_rssize_ (void)
{
  prpsinfo_t  prpsinfo;

     if ( r3_fd == 0 ) {
        sprintf(pfile, "/proc/%05d", getpid());
        r3_fd = open(pfile, O_RDONLY);
     }

     ioctl(r3_fd, PIOCPSINFO, (void *)&prpsinfo);
     return ((double)prpsinfo.pr_rssize / 64.0);
}


double r3_cpufreq_ (void)
{
   char cpufq[128];

      syssgi (SGI_GETNVRAM, "cpufreq", cpufq);
      return ((double)atoi(cpufq)*1.e9);
}


long long r3_flops_ (void)
{
   hwperf_cntr_t cnts;

      if ( r3_perf == 0 ) r3_perf = r3_flops_init_();

      ioctl(r3_fd, PIOCGETEVCTRS, (void *)&cnts);
      return ((long long)cnts.hwp_evctr[21]);
}

int r3_flops_init_ (void)
{
   int i;
   hwperf_profevctrarg_t evctr;

      for (i = 0; i < 32; i++) {
         evctr.hwp_evctrargs.hwp_evctrl[i].hwperf_spec = 0;
         evctr.hwp_ovflw_freq[i] = 0;
      }

      evctr.hwp_evctrargs.hwp_evctrl[0].hwperf_creg.hwp_mode=HWPERF_CNTEN_U;
      evctr.hwp_evctrargs.hwp_evctrl[0].hwperf_creg.hwp_ie=1;
      evctr.hwp_evctrargs.hwp_evctrl[0].hwperf_creg.hwp_ev=0;
      evctr.hwp_evctrargs.hwp_evctrl[21].hwperf_creg.hwp_mode=HWPERF_CNTEN_U;
      evctr.hwp_evctrargs.hwp_evctrl[21].hwperf_creg.hwp_ie=1;
      evctr.hwp_evctrargs.hwp_evctrl[21].hwperf_creg.hwp_ev=5;
      evctr.hwp_ovflw_sig=0;

      if ( r3_fd == 0 ) {
         sprintf(pfile, "/proc/%05d", getpid());
         r3_fd = open(pfile, O_RDONLY);
      }

      ioctl(r3_fd, PIOCENEVCTRS, (void *)&evctr);
      return 1;
}
#endif

void r3_info_init_(void);
void r3_info_print_(int *mode, int *ind, char *text, int len);
void r3_info_begin_(int *ind, char *text, int len);
void r3_info_end_(int *ind);
void r3_info_summary_(void);
void r3_barrier_(int *ind, char *text, int len);

void r3_info_init(void)
{
   r3_info_init_();
}

void r3_info_print(int mode, int ind, char *text)
{
   int len=strlen(text);
   r3_info_print_(&mode, &ind, text, len);
}

void r3_info_begin(int *ind, char *text)
{
   int len=strlen(text);
   r3_info_begin_(ind, text, len);
}

void r3_info_end(int *ind)
{
   r3_info_end_(ind);
}

void r3_info_summary(void)
{
   r3_info_summary_();
}

void r3_barrier(int *ind, char *text)
{
   int len=strlen(text);
   r3_barrier_(ind, text, len);
}
