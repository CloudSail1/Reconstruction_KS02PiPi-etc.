#!/bin/bash
let debug=0
if [[ ${debug} == 1 ]]; then
  set -x
fi

##set share file diretory
path_share="/workfs2/bes/cloudsail/workarea/boss708/Analysis/JpsiToPiKKshortAlg/JpsiToPiKKshortAlg-00-00-01/share/jobOptions_JpsiToPiKKshort.txt"
##set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
output_level=6
##if want to use exist dst file list, change to 0, default 1
is_search_dst=1
##path to search dst file
path_data="/bes3fs/offline/data/708-1/jpsi/round12/mc/dst/"
##file number to searching, set 0 to include all dst
fileNumber=0
##the data file number in one job
nDataFile=50
##the user name to query
user_name="cloudsail"
##number of job limit, default is 500
job_limit=100
##the job name with the format ${myjob}_${iJob}.txt
myjob="inc_mc"
##output directory
path_output="/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma/run/inc_mc"


##change something below ONLY when you know what you are doing

#save work path
path_work=$(pwd)

#check output path
path_output=${path_output%/}
if [ ! -e ${path_output} ] ; then
   mkdir ${path_output} -p
fi

path_dst=${path_output}/dst_list.txt
#search dst
if [[ ${is_search_dst} == 1 ]]; then
  : > ${path_dst}
  #check relative path
  if [[ ${path_data} != /* ]]; then
    path_data=${path_work}/${path_data}
  fi
  if [[ ${fileNumber} == 0 ]]; then
    find "${path_data}" -name *.dst >> ${path_dst}
  else
    find "${path_data}" -name *.dst | shuf -n ${fileNumber} >> ${path_dst}
  fi
else
  if [[ -e ${path_dst} ]]; then 
    fileNumber=$(wc -l < ${path_dst})
  else
    echo "no dst_list.txt can be use."
    return
  fi
fi #[[ ${is_search_dst} == 1 ]]

let iJob=0
let n=1
filename_job=$(printf "${myjob}_%0${#fileNumber}d" "${iJob}")
#touch first one
: > ${path_output}/${filename_job}.txt
sed -e s:path_share:${path_share}: ./run.head >> ${path_output}/${filename_job}.txt
for recFile in $(cat ${path_dst})
do
  #echo $recFile
  if [[ ! -e $recFile ]]; then 
    echo "this data file $recFile doesn't exist! "
  else
    line="\"$recFile\""
    if [[ ${n}%${nDataFile} -eq 0 || ${n} -eq ${fileNumber} ]]; then
      echo ${line} >> ${path_output}/${filename_job}.txt
      sed -e s:filename_output:${filename_job}: \
        -e s:output_level:${output_level}: \
        ./run.foot >> ${path_output}/${filename_job}.txt
      cd ${path_output}
      boss.condor ${filename_job}.txt
      sleep 0
      cd ${path_work}
      njob=$(hep_q -u "${user_name}" | tail -n 1 | awk '{print $1}')
      echo ${njob}
      while [ ${njob} -gt ${job_limit} ]
      do
        sleep 0
        njob=$(hep_q -u "${user_name}" | tail -n 1 | awk '{print $1}')
        echo ${njob}
      done
      #touch next one
      let iJob++
      filename_job=$(printf "${myjob}_%0${#fileNumber}d" "${iJob}")
      : > ${path_output}/${filename_job}.txt
      sed -e s:path_share:${path_share}: ./run.head >> ${path_output}/${filename_job}.txt
    else
      echo $line, >> ${path_output}/${filename_job}.txt
    fi
    let n++
  fi
done
#remove last noe
rm ${path_output}/${filename_job}.txt

if [[ ${debug} == 1 ]]; then
  set +x
fi
