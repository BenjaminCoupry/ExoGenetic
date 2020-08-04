using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ExoGenetic
{
    class Program
    {
        static void Main(string[] args)
        {
        }
    }
    
    //TODO, shema de genes, exprimer le genome complet...
    class Genetique
    {
        uint nbTypesBases;
        uint nbPairesChromosomes;
        List<uint> TailleChromosomes;
        uint PariteChromosomes;
        uint TailleCodon;
        uint nbAcideAmine;
        uint nbPhenoptypes;
        double probaCrossOver;
        Random r;
        public CodeGenetique Croiser(CodeGenetique[] parents)
        {
            GroupeChromosome[] code = new GroupeChromosome[nbPairesChromosomes];
            for (int i=0;i<nbPairesChromosomes;i++)
            {
                Chromosome[] paire = new Chromosome[PariteChromosomes];
                for(int j=0;j<PariteChromosomes;j++)
                {
                    paire[j] = parents[j].code[i].ExtraireChromosomeAlea(ref r, probaCrossOver);
                }
                GroupeChromosome gr = new GroupeChromosome(paire);
                code[i] = gr;
            }
            return new CodeGenetique(code);
        }
    }
    class CodeGenetique
    {
        
        public GroupeChromosome[] code;

        public CodeGenetique(GroupeChromosome[] code)
        {
            this.code = code;
        }

        public CodeGenetique(uint nbPairesChromosomes, uint nbTypesBases, List<uint> TailleChromosomes, uint PariteChromosomes)
        {
            code = new GroupeChromosome[nbPairesChromosomes];
            for(int i=0;i<nbPairesChromosomes;i++)
            {
                code[i] = new GroupeChromosome(PariteChromosomes, nbTypesBases, TailleChromosomes[i]);
            }
        }
    }
    class Chromosome
    {
        public readonly uint[] bases;
        public Chromosome(List<uint> bases_)
        {
            bases = bases_.ToArray();
        }
        public Chromosome(uint nbTypesBases, uint TailleChromosome)
        {
            bases = new UInt32[TailleChromosome];
            for(int i=0;i<TailleChromosome;i++)
            {
                bases[i] = nbTypesBases-1;//Random...
            }
        }
        public List<Codon> ExtraireCodons(uint TailleCodon)
        {
            List<Codon> retour = new List<Codon>();
            uint i0 = 0;
            while(i0+TailleCodon<bases.Length)
            {
                retour.Add(new Codon(TailleCodon, i0, this));
                i0 += TailleCodon;
            }
            return retour;
        }
        public List<Proteine> Exprimer(uint nbTypesBases, uint nbAcideAmine, uint TailleCodon)
        {
            List<Proteine> Proteines = new List<Proteine>();
            List<Codon> codons = this.ExtraireCodons(TailleCodon);
            Proteine prot = new Proteine();
            for(int i=0;i<codons.Count;i++)
            {
                uint aa = codons[i].Decoder(nbTypesBases, nbAcideAmine);
                if (aa==0)
                {
                    if(!prot.Vide())
                    {
                        Proteines.Add(prot);
                        prot = new Proteine();
                    }
                }
                else
                {
                    prot.Prolonger(aa);
                }
            }
            if (!prot.Vide())
            {
                Proteines.Add(prot);
            }
            return Proteines;
        }
    }
    class GroupeChromosome
    {
        public Chromosome[] paire;

        public GroupeChromosome(Chromosome[] paire)
        {
            this.paire = paire;
        }

        public GroupeChromosome(uint parite, uint nbTypesBases, uint TailleChromosome)
        {
            paire = new Chromosome[parite];
            for(int i=0;i<parite;i++)
            {
                paire[i] = new Chromosome(nbTypesBases, TailleChromosome);
            }
        }
        public Chromosome ExtraireChromosomeAlea(ref Random r, double probaCrossOver)
        {
            List<uint> bases = new List<uint>();
            int choix = r.Next(paire.Count());
            int i = 0;
            while(i<paire[choix].bases.Count())
            {
                bases.Add(paire[choix].bases[i]);
                i++;
                if(r.NextDouble()<probaCrossOver)
                {
                    choix = r.Next(paire.Count());
                }
            }
            return new Chromosome(bases);
        }
    }
    class Codon
    {
        uint[] bases;
        public Codon(uint TailleCodon, uint i0, Chromosome chromosome)
        {
            bases = new uint[TailleCodon];
            for(uint i=i0;i<TailleCodon+i0;i++)
            {
                bases[i - i0] = chromosome.bases[i];

            }
        }
        public uint Decoder(uint nbTypesBases, uint nbAcideAmine)
        {
            int n = bases.Length;
            long somme = 0;
            for(int i=0;i<n;i++)
            {
                somme += bases[i] * (long)Math.Pow(nbTypesBases, i);
            }
            return (uint)(somme % nbAcideAmine);
        }
    }
    class Proteine
    {
        private static double pseudo_Alea(int N)
        {
            //Retourne un double pseudo aleatoire N->R
            N = N * 58900;
            N = (N << 13) ^ N;
            N = (N * (N * N * 15731 + 789221)) + 1376312589;
            return 1.0 - (N & 0x7fffffff) / 1073741824.0;
        }
        List<uint> AcidesAmines;
        public Proteine()
        {
            AcidesAmines = new List<uint>();
        }
        public void Prolonger(uint acide)
        {
            AcidesAmines.Add(acide);
        }
        public bool Vide()
        {
            return AcidesAmines.Count == 0;
        }
        public Tuple<uint,double> getPhenotype(uint nbPhenoptypes)
        {
            uint hash = (uint)this.GetHashCode();
            return new Tuple<uint, double>(hash % nbPhenoptypes, pseudo_Alea((int)hash));
        }
    }
}
