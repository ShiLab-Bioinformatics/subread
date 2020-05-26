#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "LRMconfig.h"
#include "LRMchro-event.h"
#include "LRMfile-io.h"
#include "LRMbase-index.h"
#include "LRMhelper.h"

#define JUNCTION_CONFIRM_WINDOW 14

#define ceq(c,t) ((c)[0]==(t)[0] && (c)[1]==(t)[1])
#define c2eq(ch1, ch2, tg1, tg2) ((ceq(ch1, tg1) && ceq(ch2, tg2)) || (ceq(ch1, tg2) && ceq(ch2, tg1)) )

#define LRMis_donor_chars(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T') ||\
			    ((cc)[0]=='G' && (cc)[1]=='C')) 
				

int LRMpaired_chars(char * ch1, char * ch2){
	if (c2eq(ch1, ch2, "GC", "AG") || c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC")) {
		if ( ceq(ch1, "GC") || ceq(ch1, "CT") || ceq(ch1, "GT")) return 1;
	}
	return 0;
}


int LRMscanning_events_compare(void * arr, int l, int r){
	void ** arrr = (void **) arr;
	LRMcontext_t * context = arrr[0];
	int * event_ids = arrr[1];
	LRMevent_t * body_l = context -> event_space+event_ids[l];
	LRMevent_t * body_r = context -> event_space+event_ids[r];

	if(body_l->small_side > body_r->small_side)return 1;
	if(body_l->small_side < body_r->small_side)return -1;
	
	if(body_l->large_side > body_r->large_side)return 1;
	if(body_l->large_side < body_r->large_side)return -1;

	if(body_l->event_type > body_r->event_type) return 1;
	if(body_l->event_type < body_r->event_type) return -1;
	
	if(body_l -> indel_length > body_r -> indel_length) return -1; // same length, but L is del and R is ins -- prefer del than ins
	if(body_l -> indel_length < body_r -> indel_length) return 1;
	
	return -1;
}

void LRMscanning_events_merge(void * arr,  int start, int items, int items2){
	void ** arrr = (void **) arr;
	int * records = arrr[1];

	int read_1_ptr = start, read_2_ptr = start+items, write_ptr;
	int * merged_records = malloc(sizeof(int) * (items+items2));

	for(write_ptr=0; write_ptr<items+items2; write_ptr++){
		if((read_1_ptr >= start+items)||(read_2_ptr < start+items+items2 && LRMscanning_events_compare(arr, read_1_ptr, read_2_ptr) > 0))
			merged_records[write_ptr] = records[read_2_ptr++];
		else
			merged_records[write_ptr] = records[read_1_ptr++];
	}	
	memcpy(records + start, merged_records, sizeof(int) * (items+items2));
	free(merged_records);
}

void LRMscanning_events_exchange(void * arr, int l, int r){
	void ** arrr = (void **) arr;
	int * records = arrr[1];

	int tmpi;

	tmpi = records[l];
	records[l] = records[r];
	records[r] = tmpi;
}

int LRMevents_build_entries(LRMcontext_t  * context){
	int x1,side_i;

	for(x1=0; x1 < context->event_number; x1++){
		LRMevent_t * te = context->event_space+ x1;
		for(side_i = 0; side_i <2; side_i++){
			unsigned int sidepos = side_i?te->large_side:te->small_side;
			int * entries_list = LRMHashTableGet(context -> events_realignment, NULL+sidepos);
			//LRMprintf("INSERT ENTRY : %u -> %p ; SRC: %u ~ %u\n", sidepos, entries_list, te->small_side, te->large_side);
			if(NULL == entries_list){
				entries_list = malloc(sizeof(int) * 3);
				if(!entries_list){
					LRMprintf("ERROR: NO MEMORY CAN BE ALLOCATED!\n");
					assert(0);
				}
				entries_list[0]=2;
				entries_list[1]=0;
				LRMHashTablePut(context -> events_realignment , NULL+sidepos, entries_list);
			}
			int x2 = 0, inserted = 0;
			for(x2=1; x2< 1+ min( LRMMAX_EVENTS_PER_SITE , entries_list[0] ); x2++){
				//#warning ">>>>>>>>> COMMENT NEXT LINE <<<<<<<<"
				//if( x1 + 1 == entries_list[x2] )LRMprintf("REPEATED ENTRY: %d\n");

				if(0 == entries_list[x2]){
					entries_list[x2] = x1+1;
					if( x2 < entries_list[0] )entries_list[x2+1]=0;
					inserted = 1;
					break;
				}
			}
			if((!inserted) && entries_list[0] < LRMMAX_EVENTS_PER_SITE){
				int last_x1 = entries_list[0];
				entries_list[0] = LRMMAX_EVENTS_PER_SITE;
				entries_list = realloc(entries_list, sizeof(int)*(1+LRMMAX_EVENTS_PER_SITE));
				entries_list[last_x1] = x1+1;
				entries_list[last_x1+1] = 0;

				if(te -> small_side == 457511654 ) LRMprintf("INSERT_NEW EVENT : %d AT %u\n",x1, sidepos );
				
				LRMHashTablePut(context -> events_realignment, NULL+sidepos, entries_list);
			}
		}
	}
	return 0;
}

void LRMevents_reorder_merge_next(LRMcontext_t * context, int *order_index){
	LRMevent_t *prev_event = NULL, * new_space = malloc(sizeof(LRMevent_t) * 10000);
	int x1, new_space_size = 10000, new_space_used = 0;

	for(x1=0; x1 <=context->event_number; x1++){
		LRMevent_t *this_event = NULL;
		if(x1 < context->event_number) this_event = context->event_space+order_index[x1];
		if( x1 < context->event_number && prev_event!=NULL &&
			prev_event->large_side == this_event->large_side &&
			prev_event->small_side == this_event->small_side &&
			prev_event->event_type == this_event->event_type &&
			prev_event->indel_length == this_event->indel_length){
			prev_event -> supporting_reads ++;
		}else{
			if(new_space_size -1 < new_space_used){
				new_space_size*=1.7;
				new_space = realloc(new_space, sizeof(LRMevent_t)*new_space_size);
			}
			if(prev_event) memcpy(new_space+(new_space_used++), prev_event, sizeof(LRMevent_t));

			if(this_event){
				prev_event = this_event;
				prev_event -> supporting_reads = 1;
			}
		}
	}

	free(context -> event_space);
	context -> event_space = new_space;
	context -> event_space_size = new_space_size;
	context -> event_number = new_space_used;
}

int LRMevents_reorder(LRMcontext_t * context){
	int * order_index = malloc(context -> event_number*sizeof(int));
	int x1=0;
	while(x1<context -> event_number){
		order_index[x1]=x1;
		x1++;
	}
	void * sort_arr[2];
	sort_arr [0] = context;
	sort_arr [1] = order_index;
	
	LRMmerge_sort(sort_arr, context -> event_number, LRMscanning_events_compare, LRMscanning_events_exchange, LRMscanning_events_merge);
	//basic_sort(sort_arr, context -> event_number, LRMscanning_events_compare, LRMscanning_events_exchange);
	LRMevents_reorder_merge_next(context, order_index);
	
	if(0){
		LRMprintf("Total events : %d\n", context -> event_number);
		for(x1=0; x1<context -> event_number; x1++){
			LRMevent_t * te = context->event_space+ x1;
			if(1 || te -> small_side == 457511654){

				char pos1txt[100], pos2txt[100];
				LRMpos2txt(context, te->small_side, pos1txt);
				LRMpos2txt(context, te->large_side, pos2txt);

				LRMprintf("SORTED EVENT: TYPE: %d - INS %d %s ~ %s, nsup=%d\n", te -> event_type, te -> indel_length, pos1txt, pos2txt, te->supporting_reads);
			}
		}
	}

	free(order_index);
	return 0;
}

int LRMchro_event_new(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, LRMevent_t * new_event){
	//#warning " ================ NO INDEL EVENTS ================"
	if(new_event -> event_type == LRMEVENT_TYPE_INDEL) return 0;

	LRMthread_lock(&context -> event_space_lock);
	if(context -> event_space_size < context -> event_number + 1 ){
		context -> event_space_size *= 1.7;
		context -> event_space =realloc(context -> event_space, sizeof(LRMevent_t) * context -> event_space_size);
		if(!context -> event_space )return 1;
	}
	memcpy(context -> event_space+context -> event_number, new_event, sizeof(LRMevent_t));
	context -> event_number++;
	//LRMprintf("Total events after adding : %d\n", context -> event_number);
	LRMthread_lockrelease(&context -> event_space_lock);
	
	return 0;
}

int LRMevents_search(LRMcontext_t * context, unsigned int testing_pos, int search_for_smaller_coordinate, int * event_ids_buff){
	int * event_here = LRMHashTableGet(context -> events_realignment, NULL + testing_pos);
	if(NULL == event_here)return 0;
	int x1, ret = 0;
	for(x1 = 1; x1< 1+ min(LRMMAX_EVENTS_PER_SITE, event_here[0]); x1++){
		int event_i = event_here[x1] - 1; 
		if(0 > event_i)break;
		
		//LRMprintf("TESTING EVENT AT %u [%d] : %d \n", testing_pos, x1 , event_i);
		LRMevent_t * thise = context -> event_space + event_i;
		if((thise -> large_side ==  testing_pos && !search_for_smaller_coordinate)||
		   (thise -> small_side ==  testing_pos &&  search_for_smaller_coordinate) )
		   event_ids_buff[ret++]=event_i;
	}
	return ret;
}







#define LRMDP_score(x,y) dynamic_score_buffer[ (x)*dynamic_row_width + (y) ]
#define LRMDP_move(x,y) dynamic_movement_BEFORE_buffer[ (x)*dynamic_row_width + (y) ]

#define LRMSOFTCLIPPING_WINDOW 30
#define LRMSOFTCLIPPING_MATCHED 25

int LRMsoftclipping_moves(LRMcontext_t* context, LRMthread_context_t* thread_context, LRMread_iteration_context_t * iteration_context, char * move_buff, int moves, int bases_in_read, int to_large){
	int ii;
	
	int included_read_length = 0;
	int last_M = 0x7fffffff;
	int window_end = moves - 1;
	int window_start = 0;
	int window_MX = 0, window_M=0;


	//LRMprintf("MOVES for %s =%s\n", iteration_context -> read_name, move_buff);
	for(ii = moves -1; ii>= 0; ii--){
		if( move_buff[ii] == 'M' || move_buff[ii] == 'X' ){
			window_MX++;
			if(move_buff[ii] == 'M')
				window_M++;
		}
		if(window_MX == LRMSOFTCLIPPING_WINDOW)break;
	}
	window_start = ii;

	if(window_MX ==LRMSOFTCLIPPING_WINDOW){
		for(; window_start >=0 ; window_start--){
			if(move_buff[window_start]=='M' || move_buff[window_start] =='X'){
				window_MX ++;
				if(move_buff[window_start]=='M')
					window_M++;
			}

			if(window_MX > LRMSOFTCLIPPING_WINDOW){
				while(1){
					char nch = move_buff[window_end--];
					if(nch == 'M' || nch == 'X'){
						window_MX --;
						if(nch == 'M') window_M--;
						break;
					}
				}
			}

			//LRMprintf("M=%d, W = %d - %d, windows_MX=%d in %d\n", window_M, window_start, window_end, window_MX, bases_in_read);
			if(window_M < LRMSOFTCLIPPING_MATCHED)
				break;
		}
	}

	int smallwindow_Xs = 0;
	last_M = window_end;

	for(ii = window_end; ii>=0 && ii >= window_start; ii--){
		if(move_buff[ii] == 'M')
			last_M = ii;

		if(move_buff[ii] == 'X' && window_M < LRMSOFTCLIPPING_MATCHED){
			smallwindow_Xs++;
			if(smallwindow_Xs > 1) break;
		}
	}

	if(last_M > 0){
		//LRMprintf("M=%d, last 'M' at %d, windows_MX=%d in %d\n", window_M, last_M, window_MX, bases_in_read);
		for(ii = moves -1; ii>= last_M; ii--){
			if(move_buff[ii] == 'M' || move_buff[ii] == 'X' || move_buff[ii] == 'I')
				included_read_length++;
		}
	
		//assert( (last_M - 1) >=(bases_in_read - included_read_length));
		//LRMprintf("last_M=%d, included_read_length=%d in %d\n", last_M, included_read_length, bases_in_read );
		int Ss = (bases_in_read - included_read_length);
		if(Ss < 2 || last_M < 11){
			for(ii = last_M - 1; ii > (last_M - 1) - Ss ; ii--){
				if(ii<0){
					LRMprintf("MINUS_MOVE : %s , last_M = %d,  Ss = %d\n",  iteration_context -> read_name, last_M, Ss);
					return -1;
				}
				move_buff[ii]='S';
			}
			for(; ii >= 0; ii--)
				move_buff[ii]='.';
		}else{
			for(ii = last_M - 1; ii >= 0; ii--)
				move_buff[ii]='.';
			int splen = sprintf(move_buff + last_M - 10, "%dS", Ss);
			if(to_large){
				for(ii = 0; ii < splen/2; ii++){
					int tv;
					tv = move_buff[  last_M - 10 + ii];
					move_buff[ last_M - 10 +  ii] = move_buff[ last_M - 10 + splen -1 - ii ];
					move_buff[ last_M - 10 + splen -1 - ii ] = tv;
				}
			}
			move_buff[ last_M - 10+splen ] = '.';
		}
	}
	return 0;
} 

void LRMindel_dynamic_search_debug(LRMcontext_t* context, int * dynamic_score_buffer, char * dynamic_movement_BEFORE_buffer, int dynamic_row_width, int  dynamic_rows, int *best_offset_history){
	int rr, ii;

	LRMprintf("     ");
	for(ii=0;ii<dynamic_row_width;ii++)
		LRMprintf("  % 4d  ", ii - context -> max_dynamic_indel_length);
	LRMprintf("\n");
	
	for(rr=0; rr<dynamic_rows; rr++){
		LRMprintf("%4d | %4d ", best_offset_history?best_offset_history[ rr ]:-1, rr);
		for(ii=0;ii<dynamic_row_width;ii++){
			LRMprintf("% 4d %c  ", LRMDP_score(rr,ii), LRMDP_move(rr,ii));
		}
		LRMprintf("\n");
	}
}


void LRMtest_move_buff( LRMcontext_t* context, LRMthread_context_t* thread_context, LRMread_iteration_context_t * iteration_context, char * movtxt, int moves, int rlen){
	return;
	int ii, tmpi=-1, cumlen=0, oldopt = -1;
	int rebuild_len = 0;
	for(ii = 0; ii < moves; ii++){
		int nch = movtxt[ii];
		assert(nch);
		if(nch == '.' || nch == '/') continue;
		if(nch == 'X') nch = 'M';

		if(isdigit(nch)){
			if(tmpi<0)tmpi=0;
			tmpi=tmpi*10+(nch-'0');
		}else{
			if(tmpi<0) tmpi = 1;
			if(oldopt != nch){
				if(oldopt == 'M' || oldopt == 'I' || oldopt == 'S') rebuild_len += cumlen;
				cumlen = 0;
			}
			cumlen += tmpi;
			tmpi = -1;
			oldopt = nch;
		}
	}
	if(cumlen && (oldopt == 'M' || oldopt == 'I' || oldopt == 'S'))rebuild_len += cumlen;
	if(rebuild_len != rlen){
		LRMprintf("WRONG MOVES %s (moves = %d): %d (rebuild)  != %d (subread_dist)\t%s \n", iteration_context -> read_name, moves, rebuild_len, rlen, movtxt);
	}
}

int my_debug = 0;
int LRMdynamic_in_middle(LRMcontext_t* context, LRMthread_context_t* thread_context, LRMread_iteration_context_t * iteration_context, int last_correct_base, int first_correct_base ,  unsigned int last_correct_base_on_chro, int expected_indels){
	int moves, xx1;

	char * corrected_read = iteration_context -> read_text;//+ last_correct_base;
	int high_penalty_create_gap = 0;

	if(0) {
		char postxt[100];
		LRMpos2txt(context, last_correct_base_on_chro, postxt);
		LRMprintf("Dynamic: at %s : %d - %d ; expected = %d\n", postxt, last_correct_base, first_correct_base, expected_indels);
		for(xx1 = 0 ; xx1 < first_correct_base -  last_correct_base ; xx1++ ){
			LRMprintf("%c", corrected_read[ xx1 + last_correct_base ]);
		}
		LRMprintf("\n");
		for(xx1 = 0 ; xx1 < first_correct_base -  last_correct_base - expected_indels; xx1++){
			LRMprintf("%c", LRMgvindex_get(& context -> current_base_index,  last_correct_base_on_chro + xx1));
		}
		LRMprintf("\n");
	}

	int dynamic_rows = first_correct_base - last_correct_base;
	int trying_indel_length ;

	if(abs(expected_indels) < 10) trying_indel_length = abs(expected_indels) * 2;
	else trying_indel_length = (abs(expected_indels) -10) * 3/2 + 20;

	trying_indel_length = max(trying_indel_length, 10);
	trying_indel_length = min(trying_indel_length, LRMINDEL_DYNAMIC_CHANNEL_TOLERANCE);
	
	int best_offset_history [dynamic_rows];
	int score_match = context -> dynamic_programming_score_match;
	int score_mismatch = context -> dynamic_programming_score_mismatch;
	int score_create_gap = context -> dynamic_programming_score_create_gap * (1+high_penalty_create_gap);
	int score_extend_gap = context -> dynamic_programming_score_extend_gap;
	int dynamic_row_width = 2* trying_indel_length + 1;

	
	char postxt[100];
	LRMpos2txt(context, last_correct_base_on_chro, postxt);
	//LRMprintf("DDDDD1  SIZE=%d ; last_correct = %d; first_correct = %d ; delta=%d ; pos=%u (%s)\n", dynamic_row_width * dynamic_rows, last_correct_base, first_correct_base, expected_indels, last_correct_base_on_chro, postxt);
	int * dynamic_score_buffer = (int *)thread_context -> dynamic_programming_score_buffer;
	char * dynamic_movement_BEFORE_buffer = thread_context -> dynamic_programming_movement_buffer;
	memset(dynamic_score_buffer, 0, sizeof(int) * dynamic_row_width * dynamic_rows);
	memset(dynamic_movement_BEFORE_buffer, 0, sizeof(char) * dynamic_row_width * dynamic_rows);
	char * indel_movement_buff = (char *) thread_context -> dynamic_programming_indel_movement_buf;
	int this_movement_start = thread_context -> dynamic_programming_indel_movement_start;
	
	LRMDP_score(0,  trying_indel_length  )=0;
	int read_cursor = last_correct_base, row_i, indel_i;
	unsigned int chro_cursor = last_correct_base_on_chro;

	int last_slope_offset = 0;

	if(1){
		float slope = expected_indels *1. / dynamic_rows;
		for(xx1 = 0; xx1 < dynamic_rows; xx1++)
			best_offset_history[xx1] = 1* (int)(xx1 * slope);
	}

	for(; read_cursor < first_correct_base; read_cursor++){
		row_i = read_cursor - last_correct_base;
		int slope_offset = row_i>0?best_offset_history[row_i-1]:0;
		int last_slope_delta = slope_offset - last_slope_offset;

		for(indel_i = dynamic_row_width-1 ; indel_i >=0; indel_i --){ // negative: deletion; positive: insertion
			int testing_indel = indel_i - trying_indel_length;
			if(1){
				int score_from_del = -0x7fffffff, score_from_ins = -0x7fffffff, score_from_match = -0x7fffffff;
				int is_matched_base = toupper(corrected_read[read_cursor]) == toupper(LRMgvindex_get(& context -> current_base_index, chro_cursor - slope_offset - testing_indel));

				if(row_i>0 && (indel_i-1+ last_slope_delta)>=0 && (indel_i-1+ last_slope_delta)< dynamic_row_width && LRMDP_score(row_i-1, indel_i-1 + last_slope_delta) > -0x7ffffff0)
					score_from_ins = LRMDP_score(row_i-1, indel_i-1 + last_slope_delta) + ( (LRMDP_move(row_i-1, indel_i-1+ last_slope_delta) == 'M' || LRMDP_move(row_i-1, indel_i-1 +last_slope_delta) == 'X')?score_create_gap:score_extend_gap);

				if(testing_indel < 0 || testing_indel < row_i)if(indel_i < dynamic_row_width-1 && LRMDP_score(row_i, indel_i+1) > -0x7ffffff0)
					score_from_del = LRMDP_score(row_i, indel_i+1) + ((LRMDP_move(row_i, indel_i+1) == 'M' || LRMDP_move(row_i, indel_i+1) == 'X')?score_create_gap:score_extend_gap);

				if((indel_i+ last_slope_delta)>=0 && (indel_i+ last_slope_delta)< dynamic_row_width && (row_i ==0 || LRMDP_score(row_i-1, indel_i + last_slope_delta) > -0x7ffffff0)){
  					score_from_match =(row_i > 0 ?LRMDP_score(row_i-1, indel_i + last_slope_delta): 0)+ (is_matched_base?score_match:score_mismatch);
					if(row_i == 0 && testing_indel > 0) score_from_match += score_create_gap + (testing_indel-1) * score_extend_gap;
				}
				
				int final_score = max(score_from_del, max(score_from_ins, score_from_match));
				if(testing_indel + slope_offset > 0 && row_i < slope_offset + testing_indel){
					LRMDP_score(row_i, indel_i) = score_create_gap + (score_extend_gap-1)  * (testing_indel+slope_offset) ;
					LRMDP_move(row_i, indel_i) =  'I'; 
				}else{
					LRMDP_score(row_i, indel_i) = final_score;
					if(final_score < -0x7ffffff0) LRMDP_move(row_i, indel_i) = '?';
					else LRMDP_move(row_i, indel_i) = score_from_del == final_score?'D':((score_from_ins == final_score)?'I': ( is_matched_base ?'M':'X'));
				}
			}
		}
		last_slope_offset = slope_offset;
		chro_cursor ++;
	}

	if(0) LRMindel_dynamic_search_debug(context, dynamic_score_buffer, dynamic_movement_BEFORE_buffer, dynamic_row_width, dynamic_rows, best_offset_history);

	row_i = first_correct_base - last_correct_base - 1;
	indel_i = trying_indel_length + expected_indels -(row_i >0?best_offset_history[row_i - 1]:0);

	moves = 0;
	while(row_i >= 0 && indel_i >=0 && indel_i < dynamic_row_width){
		int slope_offset =  best_offset_history[row_i-1];
		int next_slope_offset = row_i > 1?best_offset_history[row_i-2]:0;
		int last_slope_delta = slope_offset - next_slope_offset;
		//#warning "========= DO NOT ASSERT ============="
		indel_movement_buff[ this_movement_start + moves] = LRMDP_move(row_i, indel_i);
		if(indel_movement_buff[ this_movement_start + moves]=='?')LRMprintf("Assertion_Error:%s\n", iteration_context -> read_name);
		assert(indel_movement_buff[ this_movement_start +  moves]!='?');

		if(indel_movement_buff[ this_movement_start +  moves] == 'M' || indel_movement_buff[ this_movement_start +  moves] == 'X'){
			row_i--;
			indel_i += last_slope_delta;
		} else if(indel_movement_buff[ this_movement_start +  moves] == 'D')indel_i++;
		  else {
			indel_i --;
			indel_i += last_slope_delta;
			row_i--;
		}
		moves ++;

		if(row_i < 0 && indel_i < trying_indel_length)
			for(; indel_i < trying_indel_length; indel_i++) indel_movement_buff[ this_movement_start + ( moves++ ) ] ='D';

		if(moves > max( LRMDYNAMIC_MAXIMUM_GAP_LENGTH * 15, 300 ) +  context -> max_dynamic_indel_length ){
			LRMprintf("ERROR: Dynamic programming moves more than %d\n",  max( (int)(LRMDYNAMIC_MAXIMUM_GAP_LENGTH * 15), 300 ) +  context -> max_dynamic_indel_length);
			return -1;
		}
	}
	
	indel_movement_buff[ this_movement_start +  moves]=0;
	for(row_i = 0; row_i < moves/2; row_i++){
		char tmp = indel_movement_buff[ this_movement_start +  row_i];
		indel_movement_buff[ this_movement_start +  row_i] = indel_movement_buff[ this_movement_start +  moves - row_i - 1];
		indel_movement_buff[ this_movement_start +  moves - row_i - 1] = tmp;
	}

	indel_movement_buff[ this_movement_start + (moves++)]='/';
	indel_movement_buff[ this_movement_start + moves] = 0;

	LRMtest_move_buff( context,  thread_context, iteration_context , indel_movement_buff + this_movement_start, moves , dynamic_rows);
	if(0)LRMprintf("MOVES = %s\n", indel_movement_buff + this_movement_start);

	return moves;
}

int LRMdynamic_to_ends(LRMcontext_t* context, LRMthread_context_t* thread_context, LRMread_iteration_context_t * iteration_context, int last_mapped_in_read, unsigned int last_correct_base_on_chro, int search_to_3end){
	int moves = 0;
	int high_penalty_create_gap = 0;

	int last_correct_base = search_to_3end? last_mapped_in_read : 0;
	int first_correct_base = search_to_3end? iteration_context -> read_length : last_mapped_in_read;
	int head_SS = 0, tail_SS = 0;

	int this_movement_start = thread_context -> dynamic_programming_indel_movement_start;
	//LRMprintf("SEARCH_END_01 %s : start=%d, moves=%d\n", iteration_context -> read_name , this_movement_start, moves);
	char * indel_movement_buff = (char *) thread_context -> dynamic_programming_indel_movement_buf;
	int bases_in_read_no_SS = first_correct_base - last_correct_base;
	//if(bases_in_read_no_SS<3){
	//	return sprintf(indel_movement_buff + this_movement_start,"%dS/",bases_in_read_no_SS);
	//}

	if((!search_to_3end) && first_correct_base >= LRMDYNAMIC_MAXIMUM_GAP_LENGTH){
		last_correct_base = first_correct_base - LRMDYNAMIC_MAXIMUM_GAP_LENGTH;
		head_SS = last_correct_base;
		//LRMprintf("HEAD_SKIP = %d\n", head_SS);
	}
	if(search_to_3end && iteration_context -> read_length -  last_correct_base >= LRMDYNAMIC_MAXIMUM_GAP_LENGTH){
		first_correct_base = last_correct_base + LRMDYNAMIC_MAXIMUM_GAP_LENGTH;
		tail_SS = iteration_context -> read_length - first_correct_base;
	}

	char * corrected_read = iteration_context -> read_text;// + first_correct_base;
	int bases_in_read = first_correct_base - last_correct_base;
	int trying_indel_length = LRMINDEL_DYNAMIC_CHANNEL_TOLERANCE/3;
	//if( bases_in_read > 1000 ) trying_indel_length = LRMINDEL_DYNAMIC_CHANNEL_TOLERANCE/2;
	//if( bases_in_read > 10000 ) trying_indel_length = LRMINDEL_DYNAMIC_CHANNEL_TOLERANCE;
	if(bases_in_read_no_SS<2) trying_indel_length = 1;
	
	//int best_offset_history [bases_in_read];
	int best_offset_history[bases_in_read];
	int score_match = context -> dynamic_programming_score_match;
	int score_mismatch = context -> dynamic_programming_score_mismatch;
	int score_create_gap = context -> dynamic_programming_score_create_gap * (1+high_penalty_create_gap);
	int score_extend_gap = context -> dynamic_programming_score_extend_gap;
	int dynamic_row_width = 2* trying_indel_length + 1;

	int * dynamic_score_buffer = (int *)thread_context -> dynamic_programming_score_buffer;
	char * dynamic_movement_BEFORE_buffer = thread_context -> dynamic_programming_movement_buffer;
	//LRMprintf("EEEEE1 %s: Memsize=4 * %d = %d x %d\n", iteration_context -> read_name, dynamic_row_width * bases_in_read , dynamic_row_width, bases_in_read);
	memset(dynamic_score_buffer, 0, sizeof(int) * dynamic_row_width * bases_in_read);
	memset(dynamic_movement_BEFORE_buffer, 0, sizeof(char) * dynamic_row_width * bases_in_read);
	
	LRMDP_score(0,  trying_indel_length  )=0;
	unsigned int chro_cursor ;

	int last_slope_offset = 0, read_i, indel_i, previous_base_in_read = 0;

	if(0){
		char postxt[100];
		LRMpos2txt(context, last_correct_base_on_chro, postxt);
		LRMprintf("EXTEND_UNKNOWN: %s\n", postxt);
		if(!search_to_3end){
				int bb;
				bb = corrected_read[first_correct_base];
				corrected_read[first_correct_base] = 0;
				LRMprintf("READ:           %s\n", corrected_read);
				corrected_read[first_correct_base] = bb;

				LRMprintf("CHRO: ");
				for(chro_cursor = last_correct_base_on_chro - bases_in_read - 10; chro_cursor < last_correct_base_on_chro; chro_cursor ++){
					bb = LRMgvindex_get(& context -> current_base_index, chro_cursor);
					LRMprintf("%c", bb);
				}
				LRMprintf("\n\n");
		}else{
				int bb;
				bb = corrected_read[first_correct_base];
				corrected_read[first_correct_base] = 0;
				LRMprintf("READ: %s\n", corrected_read + last_correct_base);
				corrected_read[first_correct_base] = bb;

				LRMprintf("CHRO: ");
				for(chro_cursor = last_correct_base_on_chro; chro_cursor < last_correct_base_on_chro + 10 + bases_in_read; chro_cursor ++){
					bb = LRMgvindex_get(& context -> current_base_index, chro_cursor);
					LRMprintf("%c", bb);
				}
				LRMprintf("\n\n");
		}
	}

	for(read_i = 0; read_i < bases_in_read; read_i ++){
		int this_base_in_read = search_to_3end? read_i :(bases_in_read - read_i -1);
		int this_base_value_in_read = corrected_read [ last_correct_base + this_base_in_read ];

		int slope_offset = read_i>0?best_offset_history[previous_base_in_read]:0;
		int last_slope_delta = slope_offset - last_slope_offset;
		if(0 && ! search_to_3end)LRMprintf("GET READ_BASE='%c' LAST_OFF=%d  SLP_OFF=%d  LAST_DELTA=%d\n", this_base_value_in_read, last_slope_offset, slope_offset, last_slope_delta);
		
		int thisrow_max_score = -0x7fffffff, thisrow_max_indel_from_start = -0x7fffffff;

	
		for(indel_i = search_to_3end?dynamic_row_width-1 : 0; indel_i !=(search_to_3end?-1:dynamic_row_width); indel_i+=(search_to_3end?-1:1)){
			int indel_from_start = slope_offset + ( indel_i - trying_indel_length);	// if to_3end: +:ins, -:del;  if !to_3end: +:del, -:ins
			unsigned int chro_location_after_indel = last_correct_base_on_chro + ( search_to_3end ?( read_i):(- read_i -1)) - indel_from_start; // if !to_3end: "DEL" from right is "INS" from left
			int this_base_value_in_chro = LRMgvindex_get(& context -> current_base_index, chro_location_after_indel);

			int score_from_del = -0x7fffffff, score_from_ins = -0x7fffffff, score_from_match = -0x7fffffff;
			int is_matched_base = toupper(this_base_value_in_read)== toupper(this_base_value_in_chro);


			if(search_to_3end){
					if(read_i > 0 && indel_i-1 + last_slope_delta >=0 && indel_i-1 + last_slope_delta < dynamic_row_width && read_i > 0 && LRMDP_score(read_i-1, indel_i-1 + last_slope_delta) > -0x7ffffff0)
						score_from_ins = LRMDP_score(read_i - 1, indel_i-1 + last_slope_delta) + ( (LRMDP_move(read_i - 1, indel_i-1+ last_slope_delta) == 'I')?score_extend_gap:score_create_gap);

					if( indel_from_start < 0 || indel_from_start < read_i ) if(indel_i < dynamic_row_width-1 && read_i > 0 && LRMDP_score(read_i, indel_i+1) > -0x7ffffff0)
						score_from_del = LRMDP_score(read_i, indel_i+1) + ((LRMDP_move(read_i, indel_i+1) == 'D')?score_extend_gap:score_create_gap);

					if((indel_i+ last_slope_delta)>=0 && (indel_i+ last_slope_delta)< dynamic_row_width && (read_i ==0 || LRMDP_score(read_i-1, indel_i + last_slope_delta) > -0x7ffffff0))
						score_from_match =(read_i > 0 ?LRMDP_score(read_i-1, indel_i + last_slope_delta): 0)+ (is_matched_base?score_match:score_mismatch);
					//if(read_i == 2) LRMprintf("INDEL_i=%d, F_INS=%d, F_DEL=%d, F_MAT=%d\n", indel_i, score_from_ins, score_from_del, score_from_match);
			}else{
					if( indel_from_start > 0 || ( -indel_from_start < read_i ) ) if(indel_i>0 && LRMDP_score(read_i, indel_i-1) > -0x7ffffff0)
						score_from_ins = LRMDP_score(read_i, indel_i-1) + ( (LRMDP_move(read_i, indel_i-1) == 'I')?score_extend_gap:score_create_gap);

					if(indel_i < dynamic_row_width-1 && indel_i+1 + last_slope_delta >= 0 && indel_i+1 + last_slope_delta < dynamic_row_width  && read_i > 0 && LRMDP_score(read_i - 1, indel_i+1+last_slope_delta) > -0x7ffffff0)
						score_from_del = LRMDP_score(read_i - 1, indel_i+1 + last_slope_delta) + ((LRMDP_move(read_i - 1, indel_i+1+last_slope_delta) == 'D')?score_extend_gap:score_create_gap);

					if((indel_i+ last_slope_delta)>=0 && (indel_i+ last_slope_delta)< dynamic_row_width && (read_i ==0 || LRMDP_score(read_i-1, indel_i + last_slope_delta) > -0x7ffffff0))
						score_from_match =(read_i > 0 ?LRMDP_score(read_i-1, indel_i + last_slope_delta): 0)+ (is_matched_base?score_match:score_mismatch);
					//if(read_i == 4 && !search_to_3end) // LRMprintf("INDEL_i=%d, F_INS=%d, F_DEL=%d, F_MAT=%d\n", indel_i, score_from_ins, score_from_del, score_from_match);
					//	LRMprintf("SCORE AT readi=%d, indel_i=%d = %d ;; MOVE = '%c'\n", read_i, indel_i-1+ last_slope_delta, LRMDP_score(read_i, indel_i-1+ last_slope_delta), LRMDP_move(read_i, indel_i-1+ last_slope_delta));
			}

			if(read_i == 0 && indel_from_start > 0 && !search_to_3end) score_from_match += score_create_gap + (indel_from_start-1) * score_extend_gap;
			if(read_i == 0 && indel_from_start > 0 &&  search_to_3end) score_from_match += score_create_gap + (indel_from_start-1) * score_extend_gap;

			int final_score = max(score_from_del, max(score_from_ins, score_from_match));
			//if(read_i == 4 && !search_to_3end)LRMprintf(" == READ_I %d; INDEL_I %d; INDEL_START %d; FINAL %d\n", read_i, indel_i, indel_from_start, final_score);

			if(indel_from_start < 0 && read_i < - indel_from_start && !search_to_3end){
				LRMDP_score(read_i, indel_i) = score_create_gap - score_extend_gap  * indel_from_start ;
				LRMDP_move(read_i, indel_i) =  'D';
			} else if(indel_from_start > 0 && read_i < indel_from_start && search_to_3end) {
				LRMDP_score(read_i, indel_i) = score_create_gap + score_extend_gap  * (indel_from_start -1) ;
				LRMDP_move(read_i, indel_i) =  'I';
			} else if(final_score < -0x7ffffff0){
				LRMDP_move( read_i, indel_i) = '?';
			}else{ 
				if(0 && ! search_to_3end)LRMprintf("#%d %c%c%c %c %d |", indel_i, indel_from_start?' ':'>', this_base_value_in_chro, indel_from_start?' ':'<', score_from_del == final_score?'D':((score_from_ins == final_score)?'I': ( is_matched_base ?'M':'X')), final_score);
				LRMDP_score( read_i, indel_i) = final_score;
				//LRMprintf(" !! READ_I %d; INDEL_I %d; INDEL_START %d; FINAL %d\n", read_i, indel_i, indel_from_start, final_score);
				if(final_score > thisrow_max_score && abs( indel_i - dynamic_row_width/2 ) < dynamic_row_width/3){
					thisrow_max_score = final_score;
					thisrow_max_indel_from_start = indel_from_start;
				}
				LRMDP_move(read_i, indel_i) =(score_from_del == final_score?'D':((score_from_ins == final_score)?'I': ( is_matched_base ?'M':'X')));
			}
		}
		//if(thisrow_max_indel_from_start < -0x7ffffff0) LRMprintf("ERROR READ\n");
		//else LRMprintf("      READ\n");
		assert(thisrow_max_indel_from_start > -0x7ffffff0);
		best_offset_history[this_base_in_read] = thisrow_max_indel_from_start;
		if(0 && ! search_to_3end)LRMprintf("\nSET %d-th BEST_FROM_START=%d ;  SLOPE_OFFSET=%d ;  MAX_SCORE=%d\n\n", this_base_in_read, thisrow_max_indel_from_start, slope_offset, thisrow_max_score);

		last_slope_offset = slope_offset;
		previous_base_in_read = this_base_in_read;
	}
	read_i = bases_in_read - 1;

	//LRMprintf("RESTART SEARCH: best_offset_history[read_i]=%d, best_offset_history[read_i-1]=%d\n", best_offset_history[  search_to_3end? bases_in_read - 1 : 0 ] - best_offset_history[ search_to_3end? bases_in_read - 2 : 1 ] );
	if(bases_in_read_no_SS > 2)
		indel_i = best_offset_history[  search_to_3end? bases_in_read - 1 : 0 ] - best_offset_history[ search_to_3end? bases_in_read - 2 : 1 ] + trying_indel_length;
	else
		indel_i = 1;

	int error_in = 0;
	while(read_i >= 0 && indel_i >=0 && indel_i < dynamic_row_width){
		int real_base_in_read = search_to_3end? read_i:(bases_in_read - (read_i+1));
		int slope_offset;
		int next_slope_offset;

		if(search_to_3end){
			slope_offset = (real_base_in_read>0)?best_offset_history[real_base_in_read - 1]:0;
			next_slope_offset = (real_base_in_read>1)?best_offset_history[real_base_in_read - 2]:0;
		}else{
			slope_offset = (real_base_in_read<bases_in_read - 1)?best_offset_history[real_base_in_read + 1]:0;
			next_slope_offset = (real_base_in_read<bases_in_read - 2)?best_offset_history[real_base_in_read + 2]:0;
		}

		int last_slope_delta = next_slope_offset - slope_offset;
		indel_movement_buff[ this_movement_start +  moves] = LRMDP_move(read_i, indel_i);

		//LRMprintf("%s R %d , SCORE+=%d (INDEL_i=%d)  MOVE='%c'  LAST_DELTA=%d\n", iteration_context -> read_name, read_i, LRMDP_score(read_i, indel_i), indel_i, LRMDP_move(read_i, indel_i), last_slope_delta);
		//LRMprintf("SEARCH_END_01 %s : start=%d, moves=%d\n", iteration_context -> read_name , this_movement_start, moves);
		if(indel_movement_buff[ this_movement_start +  moves]=='?'){
			error_in = 1;
			//LRMprintf("WRONG_READ_ERROR_MOVE : %s\n", iteration_context -> read_name);
			break;
		}

		if(indel_movement_buff[ this_movement_start +  moves] == 'M' || indel_movement_buff[ this_movement_start +  moves] == 'X'){
			read_i--;
			if(read_i>=0)indel_i -= last_slope_delta;
		}else {

			if(search_to_3end){
				if(indel_movement_buff[ this_movement_start +  moves] == 'D')indel_i++;
				else {
					indel_i --;
					if(read_i>=0)indel_i -= last_slope_delta;
					read_i--;
				}
			}else{
				if(indel_movement_buff[ this_movement_start +  moves] == 'I')indel_i--;
				else {
					indel_i ++;
					if(read_i>=0)indel_i -= last_slope_delta;
					read_i--;
				}
			}
		}
		moves ++;

		//if(read_i<0) LRMprintf("END_UP: indel_i = %d > %d\n", indel_i, trying_indel_length);
		if(search_to_3end){
			if(read_i < 0 && indel_i < trying_indel_length)
				for(; indel_i <trying_indel_length; indel_i++) indel_movement_buff[ this_movement_start +  (moves++)] = 'D';
		}else if(read_i < 0 && indel_i > trying_indel_length)
			for(; indel_i > trying_indel_length; indel_i--) indel_movement_buff[  this_movement_start +  (moves++)] = 'I';

		if(moves > max( LRMDYNAMIC_MAXIMUM_GAP_LENGTH * 15, 300 ) +  context -> max_dynamic_indel_length ){
			LRMprintf("ERROR: Dynamic programming moves more than %d\n",  max( (int)(LRMDYNAMIC_MAXIMUM_GAP_LENGTH * 15), 300 ) +  context -> max_dynamic_indel_length);
			return -1;
		}
	}

	//LRMprintf("SEARCH_END_02 %s : start=%d, moves=%d\n", iteration_context -> read_name , this_movement_start, moves);
	if(0){
		indel_movement_buff[this_movement_start + moves] = 0;
		LRMprintf("MOVE_ENDS [ %s E=%d] = %s\n",iteration_context  -> read_name , error_in, indel_movement_buff);
	}
	if(error_in){
			moves = sprintf(indel_movement_buff + this_movement_start, "%dS", bases_in_read);
			return moves;
	}

	indel_movement_buff[ this_movement_start +  moves]=0;
	if(!search_to_3end){
		for(read_i = 0; read_i < moves; read_i++){
			char tmp = indel_movement_buff[ this_movement_start + read_i];
			if(tmp == 'I') indel_movement_buff[ this_movement_start + read_i]='D';
			else if(tmp == 'D') indel_movement_buff[ this_movement_start + read_i]='I';
		}
	}
	
	error_in = LRMsoftclipping_moves(context,thread_context, iteration_context ,indel_movement_buff + this_movement_start, moves, bases_in_read, search_to_3end);
	if(error_in){
		LRMprintf("SOFT CLIPPING ERROR : %d of %s\n", error_in, iteration_context  -> read_name);
		return -1;
	}

	if(search_to_3end){
		for(read_i = 0; read_i < moves/2; read_i++){
			char tmp = indel_movement_buff[this_movement_start + read_i];
			indel_movement_buff[ this_movement_start + read_i] = indel_movement_buff[this_movement_start +  moves - read_i - 1];
			indel_movement_buff[ this_movement_start +  moves - read_i - 1] = tmp;
		}
	}

	if(head_SS>0){
		int iii;
		for(iii = moves-1; iii >=0; iii--){
			indel_movement_buff[ this_movement_start + iii + 10] = indel_movement_buff[ this_movement_start + iii];
		}
		for(iii = 0; iii < 10; iii++) indel_movement_buff[ this_movement_start + iii] = '.';
		moves += 10;
		iii = sprintf(indel_movement_buff + this_movement_start, "%dS",head_SS);
		indel_movement_buff[this_movement_start + iii] = '.';
	}
	if(tail_SS>0)moves += sprintf(indel_movement_buff + this_movement_start + moves, "%dS",tail_SS);
	indel_movement_buff[ this_movement_start + (moves++)]='/';
	indel_movement_buff[this_movement_start + moves] = 0;
	LRMtest_move_buff( context,  thread_context, iteration_context , indel_movement_buff + this_movement_start, moves , bases_in_read_no_SS );
	if(0)LRMprintf("MOVE_ENDS [ %s E=%d] = %s; rlen = %d\n",iteration_context  -> read_name , error_in, indel_movement_buff + this_movement_start, bases_in_read);
	return moves;
}

