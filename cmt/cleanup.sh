# echo "cleanup JpsiToPiKKshortAlg JpsiToPiKKshortAlg-00-00-01 in /workfs2/bes/cloudsail/workarea/boss708/Analysis"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtJpsiToPiKKshortAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtJpsiToPiKKshortAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=JpsiToPiKKshortAlg -version=JpsiToPiKKshortAlg-00-00-01 -path=/workfs2/bes/cloudsail/workarea/boss708/Analysis  $* >${cmtJpsiToPiKKshortAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=JpsiToPiKKshortAlg -version=JpsiToPiKKshortAlg-00-00-01 -path=/workfs2/bes/cloudsail/workarea/boss708/Analysis  $* >${cmtJpsiToPiKKshortAlgtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtJpsiToPiKKshortAlgtempfile}
  unset cmtJpsiToPiKKshortAlgtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtJpsiToPiKKshortAlgtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtJpsiToPiKKshortAlgtempfile}
unset cmtJpsiToPiKKshortAlgtempfile
return $cmtcleanupstatus

