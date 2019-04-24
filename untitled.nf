process samtools {

        tag "post-alignment processing"
        publishDir 'samtools_o'

        input:

        file samtool_bam_in from IPxi_bam

        output:

        set file('*.sorted.bam'), file('*.sorted.bam.bai'), file('*.stats.txt'), file('IPxi_reads_no.txt') into ch_samStats
        
        script:
        """
        samtools sort -@ 5 '$samtool_bam_in' -o ${samtool_bam_in.baseName}.sorted.bam
        samtools index ${samtool_bam_in.baseName}.sorted.bam ${samtool_bam_in.baseName}.sorted.bam.bai
        samtools stats ${samtool_bam_in.baseName}.sorted.bam > ${samtool_bam_in.baseName}.stats.txt
      
        cat IPxi.stats.txt | grep 'reads mapped:' | cut -f 3  > IPxi_reads_no.txt

        """

}


process bedgraph {

        tag "bedgraph"
        publishDir 'bedgraph_o'


        input:
        set file(sorted_bam_in), file(bam_index), file(stats), file(reads) from ch_samStats
        
        output:
        file '*.bedgraph' into bedgraph_out

        
        script: 
        RPM_no =  new File(${reads}).text
        """

        genomeCoverageBed -ibam "$sorted_bam_in" -bg -scale "$RPM_no" > ${sorted_bam_in.baseName}.rpm.bedgraph
        
        """

}
