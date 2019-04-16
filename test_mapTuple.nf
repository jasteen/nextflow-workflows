

ch_files = Channel.fromPath('./*.tsv')
ch_files2 = ch_files.map { tsv -> [baseName(), tsv]}
process touchFiles {

  input:
    set sample, file(tsv) from ch_files
  output:
    set sample, file('*.touched.tsv') into ch_touched

    script:

    """
    mv $tsv ${sample}.touched.tsv
    """
}

ch_temp = ch_touched.map { sample, file -> [sample, 'A', 'B', file] }


ch_temp2 = ch_temp.map{ sample, tbam, nbam, segment -> [sample, tbam, nbam, segment] }.groupTuple(by: [0,1,2])



process catfiles {
  publishDir path: './', mode: 'copy'
  input:
    set sample, tbam, nbam, file(tsv) from ch_temp2

  output:
    set sample, tbam, nbam, file("${sample}.collated.vardict.tsv") into ch_rawVardict

  script:

  myfiles = tsv.collect().join(' ')
  """
  cat ${myfiles} > ${sample}.collated.vardict.tsv
  """


}
