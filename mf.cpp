#include "BioSeq.hpp"

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::cout << "\nUsage: " << argv[0] << 
            " /path/to/data/dir max_number_mismatches min_number_matches\n" << std::endl;
        return 1;
    }

    std::fstream log_file, motif_file;
    OpenFile("MotifFinder.log", 0, log_file);
    WriteSomethingToFile(log_file, TimeStamp()+"..."+argv[0]+" "+argv[1]+" "+argv[2]+" "+argv[3]+"\n");
    BioSeq bs;
    bs.SetCountCutoff(atoi(argv[3]));
    bs.SetDelim("XXX");
    bs.SetSeqLength(50);
    bs.SetSeqAStartLength(0, 15);
    bs.SetSeqBStartLength(35, 15);
    bs.SetSeqAMisMatch(atoi(argv[2]));
    bs.SetSeqBMisMatch(atoi(argv[2]));
    bs.ReadSequences(argv[1]);
    WriteSomethingToFile(log_file, TimeStamp()+"...number of files read "+std::to_string(bs.GetNFilesRead())+"\n");
    bs.PopulateBSeqV();
    WriteSomethingToFile(log_file, TimeStamp()+"...number of sequences of length "+
            std::to_string(bs.GetSeqLength())+" found "+std::to_string(bs.GetNSeqFound())+"\n");
    bs.ClearRawSeq();
    bs.ValidateBSeqV();
    WriteSomethingToFile(log_file, TimeStamp()+"...starting search of front sequences.\n");
    bs.SearchFront(log_file);
    bs.SearchFull(log_file);
    OpenFile("Motifs.output", 0, motif_file);
    bs.WriteMotifsToFile(motif_file);
    CloseFile(motif_file);
    OpenFile("Motifs.output.detail", 0, motif_file);
    bs.WriteMotifsDetailToFile(motif_file);
    CloseFile(motif_file);
    WriteSomethingToFile(log_file, TimeStamp()+
            "...see Motifs.output and Motifs.output.detail for a list of motifs observed in no fewer than "+
            std::to_string(bs.GetCountCutoff())+" sequences.\n");
    CloseFile(log_file);
}
