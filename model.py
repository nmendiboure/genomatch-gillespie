gillespie_model = """
model rad51_recombination()
      // Parameters
      N       = 200;        // Number of free binding sites
      f       = 0.002;      // Probability of homologous association
      kon     = 0.0001;      // Association rate constant
      koff1   = 0.001;       // Base dissociation rate for an 8-nt complex
      kext1   = 0.1;        // Triplet extension probability for homologous complexes
      kext2   = 0.01;      // Triplet extension probability for heterologous complexes
      kdloop  = 0.1;        // D-loop formation rate constant
      koff2   = 0.01;      // D-loop dissociation rate constant
      kre     = 0.01;      // Recombination rate constant
    
      // Species
      S    = N;
      
      // Homologous complexes
      HM8  = 0;
      HM9  = 0;
      HM12 = 0;
      HM15 = 0;
      HM24 = 0;
      HM48 = 0;
      HM96 = 0;
      HM192 = 0;
      HM384 = 0;
      
      // Heterologous complexes
      HT8  = 0;
      HT9  = 0;
      HT12 = 0;
      HT15 = 0;
      HT24 = 0;
      HT48 = 0;
      HT96 = 0;
      HT192 = 0;
      HT384 = 0;
      
      // D-loop species
      DHM  = 0;
      DHT  = 0;
      
      // Recombined state
      R    = 0;
    
      // Association reactions
      R1: S -> HM8;  kon * f * S;
      R2: S -> HT8;  kon * (1 - f) * S;
    
      // Dissociation reactions
      R3:  HM8  -> S;  koff1 * HM8;
      R4:  HM9  -> S;  koff1 / (1.4^1) * HM9;       
      R5:  HM12 -> S;  koff1 / (1.4^2) * HM12;
      R6:  HM15 -> S;  koff1 / (1.4^3) * HM15;
      R7:  HM24 -> S;  koff1 / (1.4^3) * HM24;
      R8:  HM48 -> S;  koff1 / (1.4^3) * HM48;
      R9:  HM96 -> S;  koff1 / (1.4^3) * HM96;
      R10: HM192 -> S; koff1 / (1.4^3) * HM192;
      R11: HM384 -> S; koff1 / (1.4^3) * HM384;
      
      R12:  HT8  -> S;  koff1 * HT8;
      R13:  HT9  -> S;  koff1 / (1.4^1) * HT9;       
      R14:  HT12 -> S;  koff1 / (1.4^2) * HT12;
      R15:  HT15 -> S;  koff1 / (1.4^3) * HT15;
      R16:  HT24 -> S;  koff1 / (1.4^3) * HT24;
      R17:  HT48 -> S;  koff1 / (1.4^3) * HT48;
      R18:  HT96 -> S;  koff1 / (1.4^3) * HT96;
      R19:  HT192 -> S; koff1 / (1.4^3) * HT192;
      R20:  HT384 -> S; koff1 / (1.4^3) * HT384;
    
      // Extension reactions
      R21: HM8  -> HM9;   kext1 * HM8;
      R22: HM9  -> HM12;  kext1 * HM9;
      R23: HM12 -> HM15;  kext1 * HM12;
      R24: HM15 -> HM24;  kext1 * HM15;
      R25: HM24 -> HM48;  kext1 * HM24;
      R26: HM48 -> HM96;  kext1 * HM48;
      R27: HM96 -> HM192; kext1 * HM96;
      R28: HM192 -> HM384; kext1 * HM192;
      
      R29: HT8  -> HT9;   kext2 * HT8;
      R30: HT9  -> HT12;  kext2^2 * HT9;
      R31: HT12 -> HT15;  kext2^3 * HT12;
      R32: HT15 -> HT24;  kext2^3 * HT15;
      R33: HT24 -> HT48;  kext2^4 * HT24;
      R34: HT48 -> HT96;  kext2^4 * HT48;
      R35: HT96 -> HT192; kext2^5 * HT96;  
      R36: HT192 -> HT384; kext2^5 * HT192;
      
      // D-loop formation reactions
      R37: HM15  -> DHM;  kdloop^6 * HM15;
      R38: HM24  -> DHM;  kdloop^5 * HM24;
      R39: HM48  -> DHM;  kdloop^4 * HM48;
      R40: HM96  -> DHM;  kdloop^3 * HM96;
      R41: HM192 -> DHM;  kdloop^2 * HM192;
      R42: HM384 -> DHM;  kdloop^1 * HM384;
      
      R43: HT15  -> DHT;  kdloop^6 * HT15;
      R44: HT24  -> DHT;  kdloop^5 * HT24;
      R45: HT48  -> DHT;  kdloop^4 * HT48;
      R46: HT96  -> DHT;  kdloop^3 * HT96;
      R47: HT192 -> DHT;  kdloop^2 * HT192;
      R48: HT384 -> DHT;  kdloop^1 * HT384;
        
      // D-loop fate reactions
      R49: DHM -> S;  koff2 * DHM;
      R50: DHT -> S;  koff2 * DHT;
      R51: DHM -> R;  kre * DHM;
      R52: DHT -> R;  kre * DHT;
end

"""