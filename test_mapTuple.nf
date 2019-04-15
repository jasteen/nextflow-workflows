

ch_files = Channel.fromPath('./*.tsv')

ch_temp = ch_files.map { file -> ['S1', 'S1_FFPE.consensus.aligned.bam', 'S1_WB.consensus.aligned.bam', file] }


ch_temp2 = ch_temp.map{ sample, tbam, nbam, segment -> [sample, tbam, nbam, segment] }.groupTuple(by: [0,1,2])



process catfiles {
  publishDir path: './', mode: 'copy'
  input:
    set sample, tbam, nbam, file(tsv) from ch_temp2

  output:
    set sample, tbam, nbam, file("${sample}.collated.vardict.tsv") into ch_rawVardict

  script:

  myfiles = tsv.collec().join(' ')
  """
  cat ${myfiles} > ${sample}.collated.vardict.tsv
  """


}
