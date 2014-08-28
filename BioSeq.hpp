#if !defined(_BIOSEQ_HPP)

#define _BIOSEQ_HPP

#include "aux.hpp"
#include "dirent.h"

struct BSeq
{
    unsigned int index;   //which sequence it came from
    std::string a;
    std::string b;
};
struct Motif
{
    BSeq bs;
    std::vector < unsigned int > aQ; //measures contribution of each nuc in the motif
    std::vector < unsigned int > bQ;
    unsigned int score; // number of sequences matched
    std::vector < unsigned int > a_mis;
    std::vector < unsigned int > b_mis;
    unsigned int G_count; //to keep track of number of nuc occurances while matching
    unsigned int A_count;
};
class BioSeq
{
public:
    BioSeq(void);
    ~BioSeq(void);

    void SetDelim(const std::string &s){m_delim=s;} 
    void SetSeqLength(const unsigned int x){m_seq_length=x;}
    void SetSeqAStartLength(const unsigned int x, const unsigned int l){m_seqA_start=x;m_seqA_length=l;}
    void SetSeqAMisMatch(const unsigned int n){m_seqA_mismatches=n;}
    void SetSeqBStartLength(const unsigned int x, const unsigned int l){m_seqB_start=x;m_seqB_length=l;}
    void SetSeqBMisMatch(const unsigned int n){m_seqB_mismatches=n;}
    void SetCountCutoff(const unsigned int x){m_count_cutoff=x;}
    unsigned int GetNFilesRead(void){return m_N_files_read;}
    unsigned int GetNSeqFound(void){return m_N_BSeq_found;}
    unsigned int GetSeqLength(void){return m_seq_length;}
    unsigned int GetSeqALength(void){return m_seqA_length;}
    unsigned int GetSeqBLength(void){return m_seqB_length;}
    unsigned int GetSeqAMissMatch(void){return m_seqA_mismatches;}
    unsigned int GetSeqBMissMatch(void){return m_seqB_mismatches;}
    unsigned int GetCount(void){return m_count;}
    unsigned int GetCountCutoff(void){return m_count_cutoff;}
    Motif GetMotif(void){return m_motif;}
    std::vector < std::vector < std::string > > GetRawSequences(void){return m_raw_seq;}

    int ReadSequences(const char *path);
    void PopulateBSeqV(void);
    void ClearRawSeq(void){std::vector < std::vector < std::string > >().swap(m_raw_seq);}
    void ValidateBSeqV(void);
    bool IsNew(const std::string &str_a);
    void AddToSet(const std::string &str_a);
    bool IsMotifFront(const BSeq &bs_0);
    bool IsMotif(void);
    void SearchFront(std::fstream &myfile);
    void SearchFull(std::fstream &myfile);
    void AccumMatchtingStats(const BSeq &bs);
    void WriteMotifsToFile(std::fstream &myfile);
    void WriteMotifsDetailToFile(std::fstream &myfile);
    bool NucleotidesEqualAG(const char &n1, const char &n2);
    unsigned int Distance(const std::string &s1, const std::string &s2, const unsigned int StopIfReached);
    void AGtoAG(Motif &m);
    void AGtoO(std::string &mystr);


private:
    //input description
    std::string m_delim;           //exon delim
    unsigned int m_seq_length;
    unsigned int m_seqA_start;
    unsigned int m_seqA_length;
    unsigned int m_seqA_mismatches;  //max num of mismatches 
    unsigned int m_seqB_start;
    unsigned int m_seqB_length;
    unsigned int m_seqB_mismatches;
    unsigned int m_count_cutoff;   //only keep motifs matching at least this many sequences
    unsigned int m_N_files_read;
    unsigned int m_N_BSeq_found;

    std::vector < std::vector < std::string > > m_raw_seq;
    std::vector < std::vector < BSeq > > m_bSeqV; //main container for reshaped sequences
    Motif m_motif;
    unsigned int m_count;
    std::vector < Motif > m_motifs_a;  //string a only
    std::vector < Motif > m_motifs;    //full motifs
    std::set < std::string > m_unique_motifs; //string a only, for checking uniquenes
};

#endif
