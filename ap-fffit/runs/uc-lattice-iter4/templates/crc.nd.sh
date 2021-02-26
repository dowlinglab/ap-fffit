{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash
#$ -N {{ id }}
#$ -pe smp {{ np_global }}
#$ -r n
#$ -m ae
#$ -q long@@dowling
#$ -M bbefort@nd.edu

module purge
module load ompi
module load intel
export KMP_AFFINITY=verbose,none
export KMP_BLOCKTIME=0
export PATH=/afs/crc.nd.edu/group/maginn/group_members/Garrett_Tow/:${PATH}

{% block tasks %}
{% endblock %}
{% endblock %}

