# echo "cleanup JpsiToPiKKshortAlg JpsiToPiKKshortAlg-00-00-01 in /workfs2/bes/cloudsail/workarea/boss708/Analysis"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtJpsiToPiKKshortAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtJpsiToPiKKshortAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=JpsiToPiKKshortAlg -version=JpsiToPiKKshortAlg-00-00-01 -path=/workfs2/bes/cloudsail/workarea/boss708/Analysis  $* >${cmtJpsiToPiKKshortAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=JpsiToPiKKshortAlg -version=JpsiToPiKKshortAlg-00-00-01 -path=/workfs2/bes/cloudsail/workarea/boss708/Analysis  $* >${cmtJpsiToPiKKshortAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtJpsiToPiKKshortAlgtempfile}
  unset cmtJpsiToPiKKshortAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtJpsiToPiKKshortAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtJpsiToPiKKshortAlgtempfile}
unset cmtJpsiToPiKKshortAlgtempfile
exit $cmtcleanupstatus

