#pragma once
class Chirality
{

	int n; // hamada n parameter
	int m; //hamada m parameter

public:
	Chirality();
	Chirality(int n_new, int m_new) : n(n_new), m(m_new) {};
	~Chirality();
	int getn();
	int getm();
	void setn(int n_new);
	void setm(int m_new);
};

