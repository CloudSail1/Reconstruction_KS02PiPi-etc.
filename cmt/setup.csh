# echo "setup JpsiToPiKKshortAlg JpsiToPiKKshortAlg-00-00-01 in /workfs2/bes/cloudsail/workarea/boss708/Analysis"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtJpsiToPiKKshortAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtJpsiToPiKKshortAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=JpsiToPiKKshortAlg -version=JpsiToPiKKshortAlg-00-00-01 -path=/workfs2/bes/cloudsail/workarea/boss708/Analysis  -no_cleanup $* >${cmtJpsiToPiKKshortAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=JpsiToPiKKshortAlg -version=JpsiToPiKKshortAlg-00-00-01 -path=/workfs2/bes/cloudsail/workarea/boss708/Analysis  -no_cleanup $* >${cmtJpsiToPiKKshortAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtJpsiToPiKKshortAlgtempfile}
  unset cmtJpsiToPiKKshortAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtJpsiToPiKKshortAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtJpsiToPiKKshortAlgtempfile}
unset cmtJpsiToPiKKshortAlgtempfile
exit $cmtsetupstatus

