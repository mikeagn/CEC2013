/******************************************************************************
 * Version: 1.0.1
 * Last modified on: 27	October, 2016 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
#ifndef __SORT_IND_H__
#define __SORT_IND_H__
/** Sort indexes based on a vector of values 
 * based on: http://ideone.com/HOnvI 
 **/
#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>

enum {ASCEND, DESCEND};

template<typename T>
bool ascend_sort(std::pair<size_t, T> i, std::pair<size_t, T> j) { 
	return j.second < i.second; 
}

template<typename T>
bool descend_sort(std::pair<size_t, T> i, std::pair<size_t, T> j) {
	return i.second > j.second;
}

template<typename T>
void sortIdx(std::vector<size_t>& idx, const std::vector<T>& src, const int &dir=ASCEND)
{
	/* Resize idx to the src.size() */
	idx.resize(src.size());
	/* Create a vector of pairs (tmp) */
	std::vector< std::pair<size_t, T> > tmp;
	for (size_t i=0; i<src.size(); ++i) {
		tmp.push_back(std::pair<size_t, T> (i, src[i]));
	}

	/* Sort indeces ASCEND/DESCEND */
	if (dir == ASCEND) {
		sort(tmp.begin(), tmp.end(), ascend_sort<T>);
	} else {
		sort(tmp.begin(), tmp.end(), descend_sort<T>);
	}

	/* put sorted indeces to idx vector */
	for (size_t i=0; i<src.size(); i++){
		idx[i] = (tmp[i].first);
		//std::cout << tmp[i].first << " " << tmp[i].second << std::endl;
	}
}

#endif
