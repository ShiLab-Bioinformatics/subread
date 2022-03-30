#include<stdio.h>
#include<string.h>
#include"LRMconfig.h"
#include"LRMhelper.h"

void LRMbasic_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r)){
	int i, j;
	for(i=start; i< start + items - 1; i++)
	{
		int min_j = i;
		for(j=i + 1; j< start + items; j++)
		{
			if(compare(arr, min_j, j) > 0)
				min_j = j;
		}
		if(i!=min_j)
			exchange(arr, i, min_j);
	}
}

void LRMbasic_sort(void * arr, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r)){
	LRMbasic_sort_run(arr, 0, items, compare, exchange);
}


void LRMmerge_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	if(items > 6)
	{
		int half_point = items/2;
		LRMmerge_sort_run(arr, start, half_point, compare, exchange, merge);
		LRMmerge_sort_run(arr, start + half_point, items - half_point, compare, exchange, merge);
		merge(arr, start, half_point, items - half_point);
	}
	else
	{
		LRMbasic_sort_run(arr, start, items, compare, exchange);
	}
}
void LRMmerge_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	LRMmerge_sort_run(arr, 0, arr_size, compare, exchange, merge);
}


int LRMsmith_waterman_match_score(int is_matched){
	return is_matched?1:0;
}

int LRMsmith_waterman_linear_gap_score(int gapplen){
	return -(gapplen);
}

int LRMsmith_waterman_gap_score(int gapplen){
	return -(gapplen + 1);
}

#define LRM_sw_H(X,Y) score_buff[ (X)*(m+1) + (Y)]    // the matrix has nrow=n+1 and ncol=m+1 . Data: h00, h01, h02, ..., h0m,   h10, h11, h12, ..., h1m,   h20, h21, h22, ...
#define LRM_sw_Move(X,Y) move_buff[ (X)*(m+1) + (Y)]    // the matrix has nrow=n+1 and ncol=m+1 . Data: h00, h01, h02, ..., h0m,   h10, h11, h12, ..., h1m,   h20, h21, h22, ...
#define LRC_sw_rec_moves(M) {if(actual_moves < max_moves){moves[actual_moves] = (M) ; actual_moves++;}}
#define MOVE_Hij_HpS 0
#define MOVE_Hij_negWk 1
#define MOVE_Hij_negWl 2

// A : read; B : chromosome.
int LRMsmith_waterman(char * A, int n, char * B, int m, char * moves, int max_moves, int * score_buff_provided, int * move_buff_provided, int * high_score){
	int * score_buff = score_buff_provided, score_buff_need_free = 0, actual_moves;
	int * move_buff = move_buff_provided, move_i;
	if(high_score) *high_score=0;
	if(!score_buff) {
		int buff_size = (n +1) * (m +1);
		score_buff = malloc(buff_size * sizeof(int));
		move_buff = malloc(buff_size * sizeof(int));
		
		score_buff_need_free = 1;
	}

	// Follow the pseudo code on https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
	// Note: indices of A and B are both 1-based in the pseudo code.
	// We used the non-linear gap penalty version, hence it is a little more complicated.
	int k, l;
	for(k=0; k<= n; k++) LRM_sw_H(k,0)=0;
	for(l=0; l<= m; l++) LRM_sw_H(0,l)=0;

	int i, j, highest_i=-1, highest_j=-1, all_high=-1;
	
	for(i = 1; i <= n; i++){
		for(j = 1; j <= m; j++){
			int is_matched = A[i-1] == B[j-1];
			int Hij_HpS = LRM_sw_H(i-1,j-1) + LRMsmith_waterman_match_score(is_matched);

			int Hij_negWk = -999999999, best_k = -1;
			for(k=1; k<= i; k++){
				int testval = LRM_sw_H(i-k, j) + LRMsmith_waterman_gap_score(k);
				if(testval > Hij_negWk){
					Hij_negWk = testval;
					best_k = k;
				}
			}

			int Hij_negWl = -999999999, best_l = -1;
			for(l=1; l<= j; l++){
				// M or X
				int testval = LRM_sw_H(i, j-l) + LRMsmith_waterman_gap_score(l);
				if(testval > Hij_negWl){
					Hij_negWl = testval;
					best_l = l;
				}
			}
			int move_this = MOVE_Hij_HpS, highest_score = Hij_HpS;
			if(Hij_negWk > highest_score){
				// Insertion in "A" : i.e. an insertion when A is read and B is chro.
				highest_score = Hij_negWk;
				move_this = MOVE_Hij_negWk | (best_k<<2);
			}
			if(Hij_negWl > highest_score){
				// Insertion in "B" : i.e. a deletion when A is read and B is chro.
				highest_score = Hij_negWl;
				move_this = MOVE_Hij_negWl | (best_l <<2);
			}
			if(highest_score>0){
				LRM_sw_Move(i,j) = move_this;
				LRM_sw_H(i,j) = highest_score;
				if(highest_score > all_high) {
					all_high = highest_score;
					highest_i = i;
					highest_j = j;
				}
			}else LRM_sw_H(i,j) = 0;
		}
	}

	//LRMprintf("DP All high = %d, I, J = %d, %d HIGH\n", all_high, highest_i, highest_j);
	actual_moves = 0;
	if(highest_j>=1 && highest_i>=1){
		i = highest_i;
		j = highest_j;

		// the [here ~ end] alignment is totally worthless.
		// but we need to add "I", "D" and "S" segments to match the full lengths of A and B
		// At this moment, "A" is at the i-th base, and "B" is at the j-th base. Namely we need to add (n-i-1) bases to A after the aligned part, and add (m-j-1) bases to B after the aligned part.
		int same_length_AB_misma = min(n-i, m-j);
		int is_D_in_cigar = (n-i < m-j);// chro is longer => chro has some based deleted.
		int indel_len = is_D_in_cigar?(m - same_length_AB_misma - highest_j):(n - same_length_AB_misma - highest_i);
		while(actual_moves < indel_len)
			LRC_sw_rec_moves(is_D_in_cigar?'D':'I');
		while(actual_moves < indel_len+same_length_AB_misma) LRC_sw_rec_moves('M');
		
		while(1){
			int this_move = LRM_sw_Move(i,j) & 0x3;
			int this_score = LRM_sw_H(i,j);
			if(0 == this_score){
				// The [start ~ here] alignment is totally worthless.
				// But we need to add "I", "D" and "S" segments to match the full lengths of A and B.
				// At this moment, "A" is at the i-th base, and "B" is at the j-th base. Namely, adding i bases before the alignment of A, and adding j bases before the alignment of B.
				// Note: both i and j are 1-based indices.
				same_length_AB_misma = min(i,j) ;
				is_D_in_cigar = i<j;
				indel_len = is_D_in_cigar?(j-i):(i-j);
				for(move_i=0; move_i < same_length_AB_misma; move_i++) LRC_sw_rec_moves('M');
				for(move_i=0; move_i < indel_len; move_i++) LRC_sw_rec_moves(is_D_in_cigar?'D':'I');
				break;
			}
			if(MOVE_Hij_HpS == this_move){
				i --;
				j --;
				LRC_sw_rec_moves('M');
			}else if(MOVE_Hij_negWk == this_move){
				int best_k = LRM_sw_Move(i,j) >> 2;
				i -= best_k;
				for(move_i=0; move_i < best_k; move_i ++) LRC_sw_rec_moves('I');
			}else if(MOVE_Hij_negWl == this_move){
				int best_l = LRM_sw_Move(i,j) >> 2;
				j -= best_l;
				for(move_i=0; move_i < best_l; move_i ++) LRC_sw_rec_moves('D');
			}
		}
	}

	if(score_buff_need_free){
		free(score_buff);
		free(move_buff);
	}
	for(move_i = 0; move_i < actual_moves/2; move_i++){
		char tt = moves[move_i];
		moves[move_i] = moves[actual_moves -1 -move_i];
		moves[actual_moves -1 -move_i] = tt;
	}

	i = j = 0; 
	for(move_i = 0; move_i < actual_moves; move_i ++){
		char move = moves[move_i];
		if(move=='I') i++;
		else if(move=='D') j++;
		else {
			moves[move_i] = A[i]==B[j]?'M':'X';
			i++;
			j++;
		}
	}
	if(high_score) *high_score=all_high;
	return actual_moves;
}

// A : read; B : chromosome.
// "linear" : use the linear penalty for gaps
int LRMsmith_waterman_linear(char * A, int n, char * B, int m, char * moves, int max_moves, int * score_buff_provided, int * move_buff_provided, int * high_score){
	int * score_buff = score_buff_provided, score_buff_need_free = 0, actual_moves;
	int * move_buff = move_buff_provided, move_i;
	if(high_score) *high_score=0;
	if(!score_buff) {
		int buff_size = (n +1) * (m +1);
		score_buff = malloc(buff_size * sizeof(int));
		move_buff = malloc(buff_size * sizeof(int));
		
		score_buff_need_free = 1;
	}

	// Follow the pseudo code on https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
	// Note: indices of A and B are both 1-based in the pseudo code.
	// We used the non-linear gap penalty version, hence it is a little more complicated.
	int k, l;
	for(k=0; k<= n; k++) LRM_sw_H(k,0)=0;
	for(l=0; l<= m; l++) LRM_sw_H(0,l)=0;

	int i, j, highest_i=-1, highest_j=-1, all_high=-1;
	
	for(i = 1; i <= n; i++){
		for(j = 1; j <= m; j++){
			int is_matched = A[i-1] == B[j-1];
			int Hij_HpS = LRM_sw_H(i-1,j-1) + LRMsmith_waterman_match_score(is_matched);

			int Hij_negWk = -999999999, best_k = -1;
			for(k=1; k<= 1; k++){
				int testval = LRM_sw_H(i-k, j) + LRMsmith_waterman_linear_gap_score(1);
				if(testval > Hij_negWk){
					Hij_negWk = testval;
					best_k = k;
				}
			}

			int Hij_negWl = -999999999, best_l = -1;
			for(l=1; l<= 1; l++){
				// M or X
				int testval = LRM_sw_H(i, j-l) + LRMsmith_waterman_linear_gap_score(1);
				if(testval > Hij_negWl){
					Hij_negWl = testval;
					best_l = l;
				}
			}
			int move_this = MOVE_Hij_HpS, highest_score = Hij_HpS;
			if(Hij_negWk > highest_score){
				// Insertion in "A" : i.e. an insertion when A is read and B is chro.
				highest_score = Hij_negWk;
				move_this = MOVE_Hij_negWk | (best_k<<2);
			}
			if(Hij_negWl > highest_score){
				// Insertion in "B" : i.e. a deletion when A is read and B is chro.
				highest_score = Hij_negWl;
				move_this = MOVE_Hij_negWl | (best_l <<2);
			}

			if(highest_score>0){
				LRM_sw_H(i,j) = highest_score;
				LRM_sw_Move(i,j) = move_this;
				if(highest_score > all_high) {
					all_high = highest_score;
					highest_i = i;
					highest_j = j;
				}
			} else LRM_sw_H(i,j) = 0; // when it is 0, no further backtrace is needed.
		}
	}

	//LRMprintf("DP All high = %d, I, J = %d, %d HIGH\n", all_high, highest_i, highest_j);
	actual_moves = 0;
	if(highest_j>=1 && highest_i>=1){
		i = highest_i;
		j = highest_j;

		// the [here ~ end] alignment is totally worthless.
		// but we need to add "I", "D" and "S" segments to match the full lengths of A and B
		// At this moment, "A" is at the i-th base, and "B" is at the j-th base. Namely we need to add (n-i-1) bases to A after the aligned part, and add (m-j-1) bases to B after the aligned part.
		int same_length_AB_misma = min(n-i, m-j);
		int is_D_in_cigar = (n-i < m-j);// chro is longer => chro has some based deleted.
		int indel_len = is_D_in_cigar?(m - same_length_AB_misma - highest_j):(n - same_length_AB_misma - highest_i);
		while(actual_moves < indel_len) LRC_sw_rec_moves(is_D_in_cigar?'D':'I');
		while(actual_moves < indel_len+same_length_AB_misma) LRC_sw_rec_moves('M');
		
		while(1){
			int this_move = LRM_sw_Move(i,j) & 0x3;
			int this_score = LRM_sw_H(i,j);
			if(0 == this_score){
				// The [start ~ here] alignment is totally worthless.
				// But we need to add "I", "D" and "S" segments to match the full lengths of A and B.
				// At this moment, "A" is at the i-th base, and "B" is at the j-th base. Namely, adding i bases before the alignment of A, and adding j bases before the alignment of B.
				// Note: both i and j are 1-based indices.
				same_length_AB_misma = min(i,j) ;
				is_D_in_cigar = i<j;
				indel_len = is_D_in_cigar?(j-i):(i-j);
				for(move_i=0; move_i < same_length_AB_misma; move_i++) LRC_sw_rec_moves('M');
				for(move_i=0; move_i < indel_len; move_i++) LRC_sw_rec_moves(is_D_in_cigar?'D':'I');
				break;
			}
			if(MOVE_Hij_HpS == this_move){
				i --;
				j --;
				LRC_sw_rec_moves('M');
			}else if(MOVE_Hij_negWk == this_move){
				int best_k = LRM_sw_Move(i,j) >> 2;
				i -= best_k;
				for(move_i=0; move_i < best_k; move_i ++) LRC_sw_rec_moves('I');
			}else if(MOVE_Hij_negWl == this_move){
				int best_l = LRM_sw_Move(i,j) >> 2;
				j -= best_l;
				for(move_i=0; move_i < best_l; move_i ++) LRC_sw_rec_moves('D');
			}
		}
	}else{
		int same_length_AB_misma = min(n, m);
		int is_D_in_cigar = n < m;// chro is longer => chro has some based deleted.
		int indel_len = is_D_in_cigar?(m - same_length_AB_misma):(n - same_length_AB_misma);
		while(actual_moves < indel_len) LRC_sw_rec_moves(is_D_in_cigar?'D':'I');
		while(actual_moves < indel_len+same_length_AB_misma) LRC_sw_rec_moves('M');
	}

	if(score_buff_need_free){
		free(score_buff);
		free(move_buff);
	}

	for(move_i = 0; move_i < actual_moves/2; move_i++){
		char tt = moves[move_i];
		moves[move_i] = moves[actual_moves -1 -move_i];
		moves[actual_moves -1 -move_i] = tt;
	}

	i = j = 0; 
	for(move_i = 0; move_i < actual_moves; move_i ++){
		char move = moves[move_i];
		if(move=='I') i++;
		else if(move=='D') j++;
		else {
			moves[move_i] = A[i]==B[j]?'M':'X';
			i++;
			j++;
		}
	}
	if(high_score) *high_score=all_high;
	return actual_moves;
}
