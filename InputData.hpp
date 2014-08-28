#if !defined(_BIOSEQ_HPP)

#define _BIOSEQ_HPP

#include "auxiliary.hpp"

struct BSeq
{
    std::str a;
    std::str b;
}
class BioSeq
{
public:
    BioSeq(void);
    ~BioSeq(void);
    int ReadSequences(void);
    bool NucleotidesEqual(const std::string &n1, const std::string &n2);
    bool BSeqEqual(const BSeq &bs1, const BSeq &bs2);
    int EvaluateCandidate(void);
    int OutputMotifs(void);
    void SetCandidate(const BSeq &bs){m_bSeq_candidate=bs;} 
    void SetCountCutoff(const unsigned int x){m_count_cutoff=x;}

private:
    BSeq m_bSeq_candidate;
    unsigned int m_count;
    unsigned int m_count_cutoff;
    std::vector < BSeq > m_bSeqV;
    std::multimap < unsigned int, BSeq > m_motifs;
};

#endif
