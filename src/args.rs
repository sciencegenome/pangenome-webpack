use clap::Parser;

#[derive(Debug, Parser)]
#[clap(version)]

pub struct PangenomeArgs {
    /// please provide the kmer to be searched for the origin
    pub reads_arg: String,
    /// please provide the genome name that needs to be used for the analysis
    pub genome_arg: String,
    /// please provide the number of the threads that need to be used for the assembly
    pub thread_arg: i32,
    /// please provide the path to be searched for the strings containing the kmer
    pub protein_arg: String,
}
