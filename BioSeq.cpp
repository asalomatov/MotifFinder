#include "BioSeq.hpp"

BioSeq::BioSeq(void)
{
}
BioSeq::~BioSeq(void)
{
}
bool BioSeq::NucleotidesEqualAG(const char &n1, const char &n2)
{
    if((n1=='A' || n1=='G') && (n2=='A' || n2=='G'))
        return true;
    else
        if(n1==n2)
            return true;
    return false;
}
unsigned int BioSeq::Distance(const std::string &s1, const std::string &s2, const unsigned int StopIfReached)
{
    if(s1.length() != s2.length())
        return s1.length() + s2.length();
    unsigned int i = 0, dist =0;
    while(i < s1.length())
    {
        if(!NucleotidesEqualAG(s1[i], s2[i]))
        {
            dist++;
            if(dist==StopIfReached)
                return dist;
        }
        i++;
    }
    return dist;
}
void BioSeq::AGtoAG(Motif &m)
{
    if(m.A_count>m.G_count)
    {
        for(auto &s : m.bs.a)
            if(s=='G') s = 'A';
        for(auto &s : m.bs.b)
            if(s=='G') s = 'A';
    }
    else
    {
        for(auto &s : m.bs.a)
            if(s=='A') s = 'G';
        for(auto &s : m.bs.b)
            if(s=='A') s = 'G';
    }
}
void BioSeq::AGtoO(std::string &mystr)
{
        for(auto &s : mystr)
            if(s=='G' || s=='A') s = 'O';
}
//bool BioSeq::BSeqEqual(const BSeq &bs1, const BSeq &bs2)
//{
//    return BStrEqual(bs1.a, bs2.a) && BStrEqual(bs1.b, bs2.b);
//}
int BioSeq::ReadSequences(const char *path)
{
    DIR *dir;
    struct dirent *ent;
    int num_files = 0;
    int seq_count = 0;
    std::string path_string;
    path_string = path;
    if(path_string.back() != '/')
        path_string.append("/");
    if ((dir = opendir(path)) != NULL)
    {
        while ((ent = readdir (dir)) != NULL)
        {
            std::string file_name_string;
            file_name_string = ent->d_name;
            if(file_name_string.length()>3 && file_name_string.substr(file_name_string.length()-4, 4) == ".txt")
            {    
                std::ifstream fin((path_string + file_name_string).c_str());
                std::string line;
                if(!fin.is_open()) 
                {
                    std::cout<<"ReadSequences: could not open input file! Exiting...\n";
                    exit(1);
                }
                std::vector < std::string > raw_seq;
                while(std::getline(fin, line))
                {		
                    unsigned int i = 0;
                    unsigned int i_prev = i;
                    std::string s;
                    for(i=line.find(m_delim,0);i<line.length();i=line.find(m_delim,i_prev))
                    {
                        seq_count++;
                        s = line.substr(i_prev, i - i_prev);
                        if(s.length() >= m_seq_length) 
                            raw_seq.push_back(s);
                        i_prev = i + m_delim.length();
                       // std::cout<<i<<std::endl;
                    }
                    seq_count++;
                    s = line.substr(i_prev, line.length() - i_prev);
                    if(s.length() >= m_seq_length) 
                        raw_seq.push_back(s);
                    //m_raw_seq.push_back(line.substr(i_prev, line.length() - i_prev));
                    //std::cout<<i<<std::endl;
                }
                fin.close();
                num_files++;
                m_raw_seq.push_back(raw_seq);
            }
        }
        closedir (dir);
    } 
    else 
    {
        /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
    }
    m_N_files_read = num_files;
    m_N_BSeq_found = seq_count;
    return 0;
}
void BioSeq::PopulateBSeqV(void)
{
    if(!m_bSeqV.empty()) m_bSeqV.clear();
    std::vector < BSeq > bSeqV;
    unsigned int i = 0;
    m_N_BSeq_found = 0;
    for(auto &rs : m_raw_seq)
    {
        for(auto &s : rs)
        {
            unsigned int seqAstart = m_seqA_start;
            unsigned int seqBstart = m_seqB_start;
            while(seqBstart+m_seqB_length<=s.length())
            {
                BSeq bs;
                bs.index = i;
                bs.a = s.substr(seqAstart, m_seqA_length);
                bs.b = s.substr(seqBstart, m_seqB_length);
                bSeqV.push_back(bs);
                seqAstart++;
                seqBstart++;
            }
        }
        m_bSeqV.push_back(bSeqV);
        i++;
        m_N_BSeq_found += bSeqV.size();
        bSeqV.clear();
    }
}
void BioSeq::ValidateBSeqV(void)
{
    unsigned int num_Bseq = 0;
    for(auto &bsv : m_bSeqV)
        for(auto &bs : bsv)
        {
            num_Bseq++;
            if(bs.a.length()!=m_seqA_length || bs.b.length()!=m_seqB_length)
            {
                std::cout<<"a and b are not all of the same length, exiting...\n";
                exit(1);
            }
        }
//    std::cout<<"SeqV length : " << m_bSeqV.size() << std::endl;
//    std::cout<<"Total number of BSeq : " << num_Bseq << std::endl;
//    std::cout<<"Lengths of fragments a and b match " << m_seqA_length << " and "<< m_seqB_length << std::endl;
}
void BioSeq::SearchFull(std::fstream &myfile)
{
    int count = 0;
    for(auto &m_a : m_motifs_a)
    {
        m_motif = m_a;
        m_motif.score = 0;
        m_motif.A_count = 0;
        m_motif.G_count = 0;
        m_motif.aQ.clear();
        m_motif.bQ.clear();
        m_motif.a_mis.clear();
        m_motif.b_mis.clear();
        if(IsMotif())
        {
//            std::cout<<"motif found, score: "<<m_motif.score<<std::endl;
//            std::cout<<"motif found: "<<m_motif.bs.a+"."+m_motif.bs.b<<" "<<m_motif.bs.index<<std::endl;
//            std::cout<<"motif found A nd G counts: "<<m_motif.A_count<< " "<<m_motif.G_count<<std::endl;
//            std::cout<<"motif found aQ size: "<<m_motif.aQ.size()<<std::endl;
//            std::cout<<"motif found aQ: \n";
//            PrintVectorToConsole(m_motif.aQ, 0, m_motif.aQ.size()-1);
//            std::cout<<"motif found bQ size: "<<m_motif.bQ.size()<<std::endl;
//            std::cout<<"motif found bQ: \n";
//            PrintVectorToConsole(m_motif.bQ, 0, m_motif.bQ.size()-1);
//            std::cout<<"motif found a_mis size: "<<m_motif.a_mis.size()<<std::endl;
//            std::cout<<"motif found a_mis: \n";
//            PrintVectorToConsole(m_motif.a_mis, 0, m_motif.a_mis.size()-1);
//            std::cout<<"motif found b_mis size: "<<m_motif.b_mis.size()<<std::endl;
//            std::cout<<"motif found b_mis: \n";
//            PrintVectorToConsole(m_motif.b_mis, 0, m_motif.b_mis.size()-1);

            WriteSomethingToFile(myfile, "motif found: "+m_motif.bs.a+"."+m_motif.bs.b+
                    " score "+ std::to_string(m_motif.score)+"\n");
            m_motifs.push_back(m_motif);
        }
        if(count%1000==0)
            WriteSomethingToFile(myfile, TimeStamp()+"...processed "+
                    std::to_string(count)+" front motifs.\n");
        count++;
    }
    std::sort(m_motifs.begin(), m_motifs.end(), [](const Motif &m1, const Motif &m2){return m1.score>m2.score;});
    WriteSomethingToFile(myfile, TimeStamp()+"...finished searching, found "+
        std::to_string(m_motifs.size())+" unique  motifs.\n");
}
void BioSeq::SearchFront(std::fstream &myfile)
{
    int count = 0;
    for(auto &bsv : m_bSeqV)
        for(auto &bs : bsv)
        {
            count++;
            if(IsMotifFront(bs))
            {
                Motif m;
                m.bs = bs;
                m_motifs_a.push_back(m);
                AddToSet(bs.a);
            }
            if(count%10000==0)
                WriteSomethingToFile(myfile, TimeStamp()+"...processed "+
                        std::to_string(count)+" sequences.\n");
        }
    WriteSomethingToFile(myfile, TimeStamp()+"...finished searching front sequences, found "+
        std::to_string(m_unique_motifs.size())+" unique  motifs.\n");
}
bool BioSeq::IsMotifFront(const BSeq &bs_0)
{
    if(!IsNew(bs_0.a)) return false;
    m_count = 0;
    for(auto &bsv : m_bSeqV)
    {
        for(auto &bs : bsv)
        {
            if(m_count_cutoff+bs.index > m_N_files_read+m_count) return false;
            if(bs.index==bs_0.index)
            {
                m_count++;
                break;
            }
            if(Distance(bs.a, bs_0.a, m_seqA_mismatches+1)>m_seqA_mismatches)
                continue;
            else
            {
                m_count++;
                break;
            }
        }
        if(m_count>=m_count_cutoff) return true;
    }
    return false;
}
bool BioSeq::IsMotif(void)
{
    for(auto &bsv : m_bSeqV)
    {
        for(auto &bs : bsv)
        {
            if(m_count_cutoff+bs.index > m_N_files_read+m_motif.score) return false;
            if(bs.index==m_motif.bs.index)
            {
                m_motif.score++;
                m_motif.a_mis.push_back(0);
                m_motif.b_mis.push_back(0);
                AccumMatchtingStats(bs);
                break;
            }
            unsigned int dist_a = Distance(bs.a, m_motif.bs.a, m_seqA_mismatches+1);
            unsigned int dist_b = Distance(bs.b, m_motif.bs.b, m_seqB_mismatches+1);
            if(dist_a>m_seqA_mismatches || dist_b>m_seqB_mismatches )
                continue;
            else
            {
                m_motif.score++;
                m_motif.a_mis.push_back(dist_a);
                m_motif.b_mis.push_back(dist_b);
                AccumMatchtingStats(bs);
                break;
            }
        }
    }
    if (m_motif.score>=m_count_cutoff)
    {
        AGtoAG(m_motif);
        return true;
    }
    else
        return false;
}
void BioSeq::AccumMatchtingStats(const BSeq &bs)
{
    if(m_motif.aQ.empty())
    {
        std::vector < unsigned int > v(m_motif.bs.a.size(), 0);
        m_motif.aQ = v;
    }
    if(m_motif.bQ.empty())
    {
        std::vector < unsigned int > v(m_motif.bs.b.size(), 0);
        m_motif.bQ = v;
    }

    for(unsigned int i=0;i<bs.a.length();i++)
    {
        if(NucleotidesEqualAG(bs.a[i], m_motif.bs.a[i]))
            m_motif.aQ.at(i)++;
        if(bs.a[i] == 'G')
            m_motif.G_count++;
        else
            if(bs.a[i] == 'A')
                m_motif.A_count++;
    }
    for(unsigned int i=0;i<bs.b.length();i++)
    {
        if(NucleotidesEqualAG(bs.b[i], m_motif.bs.b[i]))
            m_motif.bQ.at(i)++;
        if(bs.b[i] == 'G')
            m_motif.G_count++;
        else
            if(bs.b[i] == 'A')
                m_motif.A_count++;
    }
}
bool BioSeq::IsNew(const std::string &str_a)
{
    std::string s = str_a;
    AGtoO(s);
    return m_unique_motifs.find(s)==m_unique_motifs.end();
}
void BioSeq::AddToSet(const std::string &str_a)
{
    std::string s = str_a;
    AGtoO(s);
    if(m_unique_motifs.find(s)!=m_unique_motifs.end())
    {
        std::cout<<"AddToSet: this motif is not new. Exiting...\n";
        exit(1);
    }
    m_unique_motifs.insert(s);
}
void BioSeq::WriteMotifsToFile(std::fstream &myfile)
{
    for(auto &m : m_motifs)
        WriteSomethingToFile(myfile, 
                std::to_string(m.score)+"\t"+m.bs.a+"."+m.bs.b+"\n");
}
void BioSeq::WriteMotifsDetailToFile(std::fstream &myfile)
{
    for(auto &m : m_motifs)
    {
        WriteSomethingToFile(myfile, std::to_string(m.score)+" | "+m.bs.a+"."+m.bs.b+" | ");
        for(auto &it : m.aQ)
            WriteSomethingToFile(myfile, std::to_string(it)+" ");
        WriteSomethingToFile(myfile, " | ");
        for(auto &it : m.bQ)
            WriteSomethingToFile(myfile, std::to_string(it)+" ");
        WriteSomethingToFile(myfile, " | ");
        for(auto &it : m.a_mis)
            WriteSomethingToFile(myfile, std::to_string(it)+" ");
        WriteSomethingToFile(myfile, " | ");
        for(auto &it : m.b_mis)
            WriteSomethingToFile(myfile, std::to_string(it)+" ");
        WriteSomethingToFile(myfile, "\n");
    }
}
