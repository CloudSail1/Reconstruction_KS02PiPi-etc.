class Comb 
{
   
 public:
     Comb() 
     {
	for(int i = 0; i < 20; i++) m_comb[i] = -1;
     }
   
     ~Comb()
     {
	;
     }
   
     int icomb(const int n) const 
     {
	return m_comb[n];
     }
   
     void setIComb(const int n, const int value) 
     {
	m_comb[n] = value;
     }
   
 private:
     int m_comb[20];
}
;

