
ch_tumor = Channel.create()
ch_normal = Channel.create()

Channel.from( [["S1_FFPE", "whatever"], ["S1_FFPE", "whatever"]] ).view()
//ch_indexedConsensusBams.choice(ch_tumor, ch_normal){ a -> a[0] =~ /FFPE$/ ? 0 : 1 }
//ch_tumor.println()
Channel.from( 1,2,3,4,5 ).view()
