#!/usr/bin/env bash
set -e


# WT1="ENST00000452863.10:c.100del ENST00000452863.10:r.100del NC_000011.10(NM_024426.6):c.100del NC_000011.10(NM_024426.6):r.100del"
# WT1="ENST00000452863.10:c.100del NC_000011.10(NM_024426.6):c.100del"

# SDHD="ENST00000375549.8:c.100del ENST00000375549.8:r.100del NC_000011.10(NM_003002.4):c.100del NC_000011.10(NM_003002.4):r.100del"
# SDHD="ENST00000375549.8:c.100del NC_000011.10(NM_003002.4):c.100del"

if [ ! $# -eq 1 ]; then
  echo "Usage: $0 name"
  exit 1
fi

name=$1

mkdir -p ${name}


# for transcript in $(python3 crossmapper.py | grep ":c." | grep -v "^NM"); do
# for transcript in $(python3 crossmapper.py | grep ":r." | grep -v "^NM"); do
# for transcript in $(python3 crossmapper.py | grep -v "^NM"); do
for transcript in $(python3 crossmapper.py); do
  fname=${name}/${transcript}.txt

  if [ -f ${fname} ];then
    echo "Skipping ${transcript}, ${fname} already exists"
    continue
  fi
  echo "Analyzing $transcript"

  if [[ ${transcript} == NM* ]]; then
    # No protein annotations for naked NM
    gtgt analyze ${transcript} --extended > ${fname}
  else
    gtgt analyze ${transcript} --protein --extended > ${fname}
  fi
done

exit


for transcript in ${WT1}; do
  echo "############### analyzing $transcript ###############"
  gtgt export ${transcript} --protein > ${transcript}.${name}.txt
done

for transcript in ${SDHD}; do
  echo "############### analyzing $transcript ###############"
  gtgt export ${transcript} --protein > ${transcript}.${name}.txt
done
