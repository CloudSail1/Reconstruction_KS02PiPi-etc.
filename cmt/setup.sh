# echo "setup JpsiToPiKKshortAlg JpsiToPiKKshortAlg-00-00-01 in /workfs2/bes/cloudsail/workarea/boss708/Analysis"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtJpsiToPiKKshortAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtJpsiToPiKKshortAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=JpsiToPiKKshortAlg -version=JpsiToPiKKshortAlg-00-00-01 -path=/workfs2/bes/cloudsail/workarea/boss708/Analysis  -no_cleanup $* >${cmtJpsiToPiKKshortAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=JpsiToPiKKshortAlg -version=JpsiToPiKKshortAlg-00-00-01 -path=/workfs2/bes/cloudsail/workarea/boss708/Analysis  -no_cleanup $* >${cmtJpsiToPiKKshortAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtJpsiToPiKKshortAlgtempfile}
  unset cmtJpsiToPiKKshortAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtJpsiToPiKKshortAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtJpsiToPiKKshortAlgtempfile}
unset cmtJpsiToPiKKshortAlgtempfile
return $cmtsetupstatus

