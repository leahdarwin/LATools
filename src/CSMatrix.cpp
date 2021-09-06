#include "CSMatrix.h"

#define ENTRY_A		1
#define ENTRY_B		0

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

void current_utc_time(struct timespec *ts) {
	#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
		clock_serv_t cclock;
		mach_timespec_t mts;
		host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
		clock_get_time(cclock, &mts);
		mach_port_deallocate(mach_task_self(), cclock);
		ts->tv_sec = mts.tv_sec;
		ts->tv_nsec = mts.tv_nsec;
	#else
		clock_gettime(CLOCK_REALTIME, ts);
	#endif
}



CSMatrix::CSMatrix(LocatingArray *locatingArray) {

	this->locatingArray = locatingArray;

	// assign factor data variable for use with grabbing column names
	this->factorData = locatingArray->getFactorData();

	// to be used when populating CSMatrix
	CSCol *csCol;
	int sum;

	// get number of factors in locating array
	int factors = locatingArray->getFactors();

	// get level counts (how many levels for each factor)
	groupingInfo = locatingArray->getGroupingInfo();

	// get the level matrix from locating array (this is the main data)
	char **levelMatrix = locatingArray->getLevelMatrix();

	// a row for every test
	rows = locatingArray->getTests();

	// initialize the data column vector
	data = new vector<CSCol*>;

	// initialize (factor + level)-to-index map
	factorLevelMap = new int*[factors];
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		factorLevelMap[factor_i] = new int[groupingInfo[factor_i]->levels];
	}

	// initialize index-to-column map (to be populated later)
	mapping = new Mapping;
	mapping->mappedTo = 0;
	mapping->mapping = NULL;

	// sum of squares of each column for normalization
	vector <float>sumOfSquares;

	// populate compressive sensing matrix column by column
	int col_i = 0;

	// Add the INTERCEPT
	csCol = new CSCol;
	sum = 0;

	csCol->factors = 0;
	csCol->setting = new FactorSetting[0];

	// populate the intercept
	for (int row_i = 0; row_i < rows; row_i++) {
		addRow(csCol);
		csCol->dataP[row_i] = ENTRY_A;
		sum += 1;
	}

	// push into vector
	data->push_back(csCol);
	sumOfSquares.push_back(sum);
	col_i++;

	// go over all 1-way interactions
	for (char factor_i = 0; factor_i < factors; factor_i++) {

		for (char level_i = 0; level_i < groupingInfo[factor_i]->levels; level_i++) {

			addOneWayInteraction(factor_i, level_i, levelMatrix, sumOfSquares);

			// set factor level map
			factorLevelMap[factor_i][level_i] = col_i - 1; // subtract 1 for INTERCEPT

			// increment column index
			col_i++;

		}

	}

	int lastOneWay_i = col_i;

	// initialize 1st level mapping
	mapping->mapping = new Mapping*[col_i - 1];

	cout << "Adding t-way interactions" << endl;
	addTWayInteractions(csCol, col_i - 1, col_i, locatingArray->getT(),
		mapping->mapping, sumOfSquares, groupingInfo, levelMatrix);

	cout << "Went over " << col_i << " columns" << endl;

	// check coverability
	for (int col_i = 0; col_i < getCols(); col_i++) {
		csCol = data->at(col_i);
		csCol->coverable = checkColumnCoverability(csCol);

		if (!csCol->coverable) {
			cout << "Not coverable: " << getColName(csCol) << endl;
		}
	}

	// perform sqaure roots on sum of squares
	for (int col_i = 0; col_i < getCols(); col_i++) {
		sumOfSquares[col_i] = sqrt(sumOfSquares[col_i]);
	}

	// normalize all columns
	for (int col_i = 0; col_i < getCols(); col_i++) {
		csCol = data->at(col_i);
		for (int row_i = 0; row_i < getRows(); row_i++) {
		//	csCol->data[row_i] = csCol->data[row_i] / sumOfSquares[col_i];
		}
		data->at(col_i) = csCol;
	}

	cout << "Finished constructing CS Matrix" << endl;

}

void CSMatrix::addRow(CSCol *csCol) {
	csCol->dataVector.push_back(0);
	csCol->dataP = &csCol->dataVector[0];
}

void CSMatrix::remRow(CSCol *csCol) {
	csCol->dataVector.pop_back();
	csCol->dataP = &csCol->dataVector[0];
}

void CSMatrix::addTWayInteractions(CSCol *csColA, int colBMax_i, int &col_i, int t,
		Mapping **mapping, vector <float>&sumOfSquares, GroupingInfo **groupingInfo, char **levelMatrix) {

	CSCol *csCol, *csColB, *csColC;
	int sum;

	// the offset is 1 for the INTERCEPT
	int colBOffset = 1;

	int colCMax_i = 0;

	// start after INTERCEPT, iterate over all main effects until colBMax_i
	for (int colB_i = 0; colB_i < colBMax_i; colB_i++) {

		// grab pointer to column
		csColB = data->at(colB_i + colBOffset);

		// the next 2 lines move colCMax_i to the 1st column with the same factor as csColB
		csColC = data->at(colCMax_i + colBOffset);
		if (csColB->setting[0].factor_i > csColC->setting[0].factor_i) colCMax_i = colB_i;

		// set mapping for this column
		mapping[colB_i] = new Mapping;
		mapping[colB_i]->mapping = (t > 1 ? new Mapping*[colCMax_i] : NULL);
		mapping[colB_i]->mappedTo = (csColA->factors > 0 ? col_i : colB_i + colBOffset);
		Mapping *groupMapping = mapping[colB_i];

		int colsInGroupB = 1;
		char groupIndexB = -1;
		char levelIndexB = csColB->setting[0].index;

		// check if this column is grouped
		if (groupingInfo[csColB->setting[0].factor_i]->grouped) {
			groupIndexB = groupingInfo[csColB->setting[0].factor_i]->levelGroups[levelIndexB];

			for (int level_i = levelIndexB + 1;
				level_i < groupingInfo[csColB->setting[0].factor_i]->levels &&
				groupingInfo[csColB->setting[0].factor_i]->levelGroups[level_i] == groupIndexB; level_i++) {

				colsInGroupB++;
				colB_i++;

				// set to group mapping unless main effect (no grouping on main effects)
				if (csColA->factors > 0) {
					mapping[colB_i] = groupMapping;
				} else {
					mapping[colB_i] = new Mapping;
					mapping[colB_i]->mapping = groupMapping->mapping;
					mapping[colB_i]->mappedTo = colB_i + colBOffset;
				}
			}

			csColB = data->at(colB_i + colBOffset);
		}

		// create new column for CS Matrix
		csCol = new CSCol;
		for (int row_i = 0; row_i < rows; row_i++) {
			addRow(csCol);
		}

		// set the headers from the combining columns
		csCol->factors = csColA->factors + 1;
		csCol->setting = new FactorSetting[csCol->factors];

		// copy previous factors
		for (int setting_i = 0; setting_i < csColA->factors; setting_i++) {
			csCol->setting[setting_i] = csColA->setting[setting_i];
		}

		// assign new setting
		csCol->setting[csColA->factors].grouped = (groupIndexB != -1);			// is 1st factor grouped
		csCol->setting[csColA->factors].factor_i = csColB->setting[0].factor_i;	// 1st factor
		csCol->setting[csColA->factors].index = levelIndexB;					// set level of 1st factor
		csCol->setting[csColA->factors].levelsInGroup = colsInGroupB;			// set last level if group

		// populate the actual column data of CS matrix
		sum = populateColumnData(csCol, levelMatrix, 0, rows);

		// add new CS column to matrix if needed
		bool colAddedToMatrix = false;
		if (csCol->factors > 1) {
			// push into vector
			data->push_back(csCol);
			sumOfSquares.push_back(sum);

			col_i++;

			colAddedToMatrix = true;
		}

		if (t > 1) {
			// recursive call
			addTWayInteractions(csCol, colCMax_i, col_i, t - 1,
				mapping[colB_i]->mapping, sumOfSquares, groupingInfo, levelMatrix);
		}

		// deallocate the column if not added to CS matrix
		if (!colAddedToMatrix) {
			delete[] csCol->setting;
			delete csCol;
		}
	}

}

int CSMatrix::populateColumnData(CSCol *csCol, char **levelMatrix, int row_top, int row_len) {
	int sum = 0;

	// populate every row
	bool rowData; // start with true, perform AND operation
	for (int row_i = row_top; row_i < row_top + row_len; row_i++) {

		rowData = true; // start with true, perform AND operation

		for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {

			rowData &= levelMatrix[row_i][csCol->setting[setting_i].factor_i] >=
				csCol->setting[setting_i].index;

			rowData &= levelMatrix[row_i][csCol->setting[setting_i].factor_i] <
				csCol->setting[setting_i].index + csCol->setting[setting_i].levelsInGroup;

		}

		// AND operation
		csCol->dataP[row_i] = (rowData ? ENTRY_A : ENTRY_B);

		// add to sum of squares
		sum += csCol->dataP[row_i] * csCol->dataP[row_i];
	}

	return sum;
}

void CSMatrix::reorderRows(int k, int c) {
	FactorSetting *settingToResample;
	int nPaths;
	long long int score;
	int cols = getCols();

	// check advanced
	CSCol **array = new CSCol*[cols];
	for (int col_i = 0; col_i < cols; col_i++) {
		array[col_i] = data->at(col_i);
	}

	int coverableMin = sortByCoverable(array, 0, getCols() - 1);
	cout << "Coverable columns begin at: " << coverableMin << endl;

	int tWayMin = sortByTWayInteraction(array, coverableMin, getCols() - 1);
	cout << "t-way interactions begin at: " << tWayMin << endl;

	long long int *rowContributions = new long long int[rows];

	Path *path = new Path;
	path->entryA = NULL;
	path->entryB = NULL;
	path->min = coverableMin;
	path->max = getCols() - 1;

	while (true) {
		nPaths = 0;
		for (int row_i = 0; row_i < rows; row_i++) rowContributions[row_i] = 0;
		pathSort(array, path, 0, nPaths, NULL);

		score = 0;
		settingToResample = NULL;
		pathLAChecker(array, path, path, 0, k, score, settingToResample, rowContributions);
		minCountCheck(array, c, score, settingToResample, rowContributions);

		cout << "Score: " << score;

		int swaps = 0;
		while (true) {

			// find the max contribution out of place
			int outOrderRow_i = -1;
			for (int row_i = 1; row_i < rows; row_i++) {
				if (rowContributions[row_i - 1] < rowContributions[row_i]) {
					outOrderRow_i = row_i;
					break;
				}
			}

			if (outOrderRow_i != -1) {
				int row_i2 = outOrderRow_i;
				for (int row_i = row_i2 + 1; row_i < rows; row_i++) {
					if (rowContributions[row_i] >= rowContributions[row_i2]) {
						row_i2 = row_i;
					}
				}

				int row_i1 = -1;
				for (int row_i = 0; row_i < rows; row_i++) {
					if (rowContributions[row_i2] > rowContributions[row_i]) {
						row_i1 = row_i;
						break;
					}
				}

				// perform swap if necessary
				swapRows(array, row_i1, row_i2);
				long long int tempContribution = rowContributions[row_i1];
				rowContributions[row_i1] = rowContributions[row_i2];
				rowContributions[row_i2] = tempContribution;
				swaps++;
			} else {
				break;
			}

		}

		cout << "\tSwaps: " << swaps << endl;
		if (swaps == 0) break;
	}

	for (int row_i = 0; row_i < rows; row_i++) cout << row_i << "\t" << rowContributions[row_i] << endl;

	delete rowContributions;
}

void CSMatrix::reorderRows(int k, int c, VectorXf ** response) {
	FactorSetting *settingToResample;
	int nPaths;
	long long int score;
	int cols = getCols();

	float *responseData = (*response)->getData();

	cout << "CHECK RESPONSE BEFORE: " << responseData[0] << endl;

	// check advanced
	CSCol **array = new CSCol*[cols];
	for (int col_i = 0; col_i < cols; col_i++) {
		array[col_i] = data->at(col_i);
	}

	int coverableMin = sortByCoverable(array, 0, getCols() - 1);
	cout << "Coverable columns begin at: " << coverableMin << endl;

	int tWayMin = sortByTWayInteraction(array, coverableMin, getCols() - 1);
	cout << "t-way interactions begin at: " << tWayMin << endl;

	long long int *rowContributions = new long long int[rows];

	Path *path = new Path;
	path->entryA = NULL;
	path->entryB = NULL;
	path->min = coverableMin;
	path->max = getCols() - 1;

	while (true) {
		nPaths = 0;
		for (int row_i = 0; row_i < rows; row_i++) rowContributions[row_i] = 0;
		pathSort(array, path, 0, nPaths, NULL);

		score = 0;
		settingToResample = NULL;
		pathLAChecker(array, path, path, 0, k, score, settingToResample, rowContributions);
		minCountCheck(array, c, score, settingToResample, rowContributions);

		cout << "Score: " << score;

		int swaps = 0;
		while (true) {

			// find the max contribution out of place
			int outOrderRow_i = -1;
			for (int row_i = 1; row_i < rows; row_i++) {
				if (rowContributions[row_i - 1] < rowContributions[row_i]) {
					outOrderRow_i = row_i;
					break;
				}
			}

			if (outOrderRow_i != -1) {
				int row_i2 = outOrderRow_i;
				for (int row_i = row_i2 + 1; row_i < rows; row_i++) {
					if (rowContributions[row_i] >= rowContributions[row_i2]) {
						row_i2 = row_i;
					}
				}

				int row_i1 = -1;
				for (int row_i = 0; row_i < rows; row_i++) {
					if (rowContributions[row_i2] > rowContributions[row_i]) {
						row_i1 = row_i;
						break;
					}
				}

				// perform swap if necessary
				swapRows(array, row_i1, row_i2);
				long long int tempContribution = rowContributions[row_i1];
				rowContributions[row_i1] = rowContributions[row_i2];
				rowContributions[row_i2] = tempContribution;

				//swap response data
				float tempResponse = responseData[row_i1];
				(*response)->getData()[row_i1] = responseData[row_i2];
				(*response)->getData()[row_i2] = tempResponse;

				swaps++;
			} else {
				break;
			}

		}

		cout << "\tSwaps: " << swaps << endl;
		if (swaps == 0) break;
	}

	for (int row_i = 0; row_i < rows; row_i++) cout << row_i << "\t" << rowContributions[row_i] << endl;

	delete rowContributions;
}

void CSMatrix::minCountCheck(CSCol **array, int c,
		long long int &score, FactorSetting *&settingToResample, long long int *rowContributions) {
	int cols = getCols();

	int *count = new int[cols];

	for (int col_i = 0; col_i < cols; col_i++) {
		count[col_i] = 0;

		for (int row_i = 0; row_i < rows && count[col_i] < c; row_i++) {
			if (array[col_i]->coverable && array[col_i]->dataP[row_i] == ENTRY_A) {
				count[col_i]++;

				// add row contributions
				if (rowContributions != NULL) {
					rowContributions[row_i]++;
				}
			}
		}
	}

	for (int col_i = 0; col_i < cols; col_i++) {
		if (array[col_i]->coverable && count[col_i] < c) {
//			cout << "Below c requirement (" << count[col_i] << "): " << getColName(array[col_i]) << endl;
			score += c - count[col_i];
			if (settingToResample == NULL) {
				// randomly choose a setting in the column to resample
				settingToResample = &array[col_i]->setting[rand() % array[col_i]->factors];
			}
		}
	}

	delete[] count;
}

void CSMatrix::exactFix() {

	// create a work array
	CSCol **array = new CSCol*[getCols()];
	for (int col_i = 0; col_i < getCols(); col_i++) {
		array[col_i] = data->at(col_i);
	}

	long long int score = 0;

	smartSort(array, 0);
	score = getArrayScore(array);

	cout << "Original linear LA Score: " << score << endl;

	if (locatingArray->getNConGroups() == 0) {
		while (score > 0) {
			addRowFix(array, score);
		}

		cout << "Complete LA created with score: " << score << endl;
		cout << "Rows: " << getRows() << endl;
	} else {
		cout << "Unable to perform fixla because constraints were found. Remove the constraints to perform this operation!" << endl;
	}

	delete[] array;

}

void CSMatrix::performCheck(int k, int c) {

	int nConGroups = locatingArray->getNConGroups();
	ConstraintGroup **conGroups = locatingArray->getConGroups();
	long long int score = 0;
	struct timespec start;
	struct timespec finish;
	float elapsedTime;

	// create a work array
	CSCol **array = new CSCol*[getCols()];
	for (int col_i = 0; col_i < getCols(); col_i++) {
		array[col_i] = data->at(col_i);
	}

	int coverableMin = sortByCoverable(array, 0, getCols() - 1);
	cout << "Coverable columns begin at: " << coverableMin << endl;

	int tWayMin = sortByTWayInteraction(array, coverableMin, getCols() - 1);
	cout << "t-way interactions begin at: " << tWayMin << endl;

	Path *path = new Path;
	path->entryA = NULL;
	path->entryB = NULL;
	path->min = coverableMin;
	path->max = getCols() - 1;

	list <Path*>pathList;

	int nPaths = 0;

	// grab initial time
	current_utc_time( &start);
	pathSort(array, path, 0, nPaths, &pathList);
	// check current time
	current_utc_time( &finish);
	// get elapsed seconds
	elapsedTime = (finish.tv_sec - start.tv_sec);
	// add elapsed nanoseconds
	elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	cout << "Elapsed for path sort: " << elapsedTime << endl;

	cout << "Memory check: nPaths = " << nPaths << " of size " << sizeof(Path) << endl;
	cout << "Unfinished paths: " << pathList.size() << endl;

	// grab initial time
	current_utc_time( &start);
	FactorSetting *settingToResample = NULL;
	pathLAChecker(array, path, path, 0, k, score, settingToResample, NULL);
	minCountCheck(array, c, score, settingToResample, NULL);
	// check current time
	current_utc_time( &finish);
	// get elapsed seconds
	elapsedTime = (finish.tv_sec - start.tv_sec);
	// add elapsed nanoseconds
	elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	cout << "Path and min count LA Score: " << score << endl;
	cout << "Elapsed for path and min count check: " << elapsedTime << endl;

	deletePath(path);

	// grab initial time
	current_utc_time(&start);
	smartSort(array, 0);
	// check current time
	current_utc_time(&finish);
	// get elapsed seconds
	elapsedTime = (finish.tv_sec - start.tv_sec);
	// add elapsed nanoseconds
	elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	cout << "Elapsed for smart sort: " << elapsedTime << endl;

	// grab initial time
	current_utc_time(&start);
	score = getArrayScore(array);
	minCountCheck(array, c, score, settingToResample, NULL);
	// check current time
	current_utc_time(&finish);
	// get elapsed seconds
	elapsedTime = (finish.tv_sec - start.tv_sec);
	// add elapsed nanoseconds
	elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	cout << "Weird linear check LA Score (should not match other scores): " << score << endl;
	cout << "Elapsed for linear check: " << elapsedTime << endl;

	/* BRUTE FORCE */
	cout << "Performing brute force check... this could take awhile" << endl;
	// grab initial time
	current_utc_time(&start);
	score = getBruteForceArrayScore(array, k);
	minCountCheck(array, c, score, settingToResample, NULL);
	// check current time
	current_utc_time(&finish);
	// get elapsed seconds
	elapsedTime = (finish.tv_sec - start.tv_sec);
	// add elapsed nanoseconds
	elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	cout << "Brute force and min count LA Score (should match path score): " << score << endl;
	cout << "Elapsed for brute force check: " << elapsedTime << endl;

	cout << "Checking for constraint violations..." << endl;
	for (int row_i = 0; row_i < rows; row_i++) {
		for (int conGroup_i = 0; conGroup_i < nConGroups; conGroup_i++) {
			if (!conGroups[conGroup_i]->getResult(row_i)) {
				cout << "Constraint group " << conGroup_i << " violated in row " << row_i << endl;
			}
		}

	}
	cout << "Done!" << endl;

	delete[] array;

}

void CSMatrix::autoFindRows(int k, int c, int startRows) {

	int iters = 1000;
	int cols = getCols();

	int upperBound = startRows;
	int lowerBound = 1;

	// check advanced
	CSCol **array = new CSCol*[cols];
	for (int col_i = 0; col_i < cols; col_i++) {
		array[col_i] = data->at(col_i);
	}

	int twoWayMin;
	for (twoWayMin = 0; twoWayMin < cols; twoWayMin++) {
		if (array[twoWayMin]->factors > 1) break;
	}
	cout << "Two-Way Min: " << twoWayMin << endl;

	int factors = locatingArray->getFactors();
	int nPaths = 0;
	long long int score;
	FactorSetting *settingToResample = NULL;

	Path *path = new Path;
	path->entryA = NULL;
	path->entryB = NULL;
//	path->min = twoWayMin;
	path->min = 0;
	path->max = getCols() - 1;

	// add more rows to reach total count
	resizeArray(array, startRows);

	// use a binary search to find the correct value
	while (true) {

		// check if it finds a proper array once in 5 times
		bool testPassed = false;
		for (int i = 0; i < 5; i++) {

			randomizeArray(array);

			list <Path*>pathList;
			pathList.push_front(path);

			score = 0;
			randomizePaths(array, settingToResample, path, 0, k, c, score, &pathList, iters);

			cout << "Score: " << score << endl;

			if (settingToResample == NULL) {
				testPassed = true;
				break;
			} else if (score > 100) {
				// do not try 5 times because the score is greater than 100
				break;
			}
		}

		// reset upper / lower bounds
		if (testPassed) {
			upperBound = rows;
		} else {
			lowerBound = rows + 1;
		}

		// check if the upper and lower bounds match
		if (upperBound < lowerBound)
			cout << "Our bounds messed up :(" << endl;
		if (upperBound <= lowerBound) break;

		// calculate a median row count
		int medRowCount = 1 * (int)((lowerBound + upperBound) / 2.0);

		// resize the array
		resizeArray(array, medRowCount);

	}

	cout << "Finished with array containing (" << lowerBound << ":" << upperBound << ") rows!" << endl;

	deletePath(path);
	delete[] array;

}

void CSMatrix::randomFix(int k, int c, int totalRows) {

	int iters = 1000;
	int cols = getCols();

	// check advanced
	CSCol **array = new CSCol*[cols];
	for (int col_i = 0; col_i < cols; col_i++) {
		array[col_i] = data->at(col_i);
	}

	int coverableMin = sortByCoverable(array, 0, cols - 1);
	cout << "Coverable columns begin at: " << coverableMin << endl;

	int tWayMin = sortByTWayInteraction(array, coverableMin, cols - 1);
	cout << "t-way interactions begin at: " << tWayMin << endl;

	int factors = locatingArray->getFactors();
	int nPaths = 0;
	long long int score;
	FactorSetting *settingToResample = NULL;

	Path *path = new Path;
	path->entryA = NULL;
	path->entryB = NULL;
	path->min = coverableMin;
	path->max = getCols() - 1;

	// add more rows to reach total count
	resizeArray(array, totalRows);

	list <Path*>pathList;
	pathList.push_front(path);

	score = 0;
	randomizePaths(array, settingToResample, path, 0, k, c, score, &pathList, iters);

	minCountCheck(array, c, score, settingToResample, NULL);

	cout << "Score: " << score << endl;

	deletePath(path);

}

void CSMatrix::systematicRandomFix(int k, int c, int initialRows, int minChunk) {
	int chunk = initialRows;
	int finalizedRows = rows;
	int totalRows = (finalizedRows + chunk);
	int cols = getCols();

	// check advanced
	CSCol **array = new CSCol*[cols];
	for (int col_i = 0; col_i < cols; col_i++) {
		array[col_i] = data->at(col_i);
	}

	int coverableMin = sortByCoverable(array, 0, cols - 1);
	cout << "Coverable columns begin at: " << coverableMin << endl;

	int tWayMin = sortByTWayInteraction(array, coverableMin, cols - 1);
	cout << "t-way interactions begin at: " << tWayMin << endl;

	Path *path = new Path;
	path->entryA = NULL;
	path->entryB = NULL;
	path->min = coverableMin;
	path->max = getCols() - 1;

	list <Path*>pathList;

	int nPaths = 0;
	pathSort(array, path, 0, nPaths, &pathList);
	cout << "nPaths: " << nPaths << " of size " << sizeof(Path) << endl;
	cout << "Unfinished paths: " << pathList.size() << endl;

	long long int score = 0;
	struct timespec start;
	struct timespec finish;
	float elapsedTime;

	// grab initial time
	current_utc_time( &start);

	FactorSetting *settingToResample = NULL;
	pathLAChecker(array, path, path, 0, k, score, settingToResample, NULL);
	minCountCheck(array, c, score, settingToResample, NULL);

	// check current time
	current_utc_time( &finish);
	// get elapsed seconds
	elapsedTime = (finish.tv_sec - start.tv_sec);
	// add elapsed nanoseconds
	elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	cout << "Elapsed: " << elapsedTime << endl;

	cout << "Score: " << score << endl;
	int factors = locatingArray->getFactors();

	while (settingToResample != NULL) {
		resizeArray(array, totalRows);

		randomizePaths(array, settingToResample, path, finalizedRows, k, c, score, &pathList, 1000);
		finalizedRows = rows;

		chunk -= chunk / 2;
		if (chunk < minChunk) chunk = minChunk;
		totalRows += chunk;
	}

}

void CSMatrix::randomizePaths(CSCol **array, FactorSetting *&settingToResample, Path *path, int row_top, int k, int c, long long int &score, list <Path*>*pathList, int iters) {

	int cols = getCols();
	long long int newScore;
	int factor_i, nPaths;
	ConstraintGroup *conGroup;
	char **levelMatrix = locatingArray->getLevelMatrix();
	int factors = locatingArray->getFactors();

	// custom factors to resample
	int nCustomFactors = 0;
	int *customFactorIndeces = new int[nCustomFactors];
//	customFactorIndeces[0] = 16;
//	customFactorIndeces[1] = 17;

	// allocate memory for saving old levels of locating array
	char **oldLevels = new char*[rows];
	for (int row_i = 0; row_i < rows; row_i++) {
		oldLevels[row_i] = new char[factors];
	}

	// sort paths
	for (std::list<Path*>::iterator it = pathList->begin(); it != pathList->end(); it++) {
		nPaths = 0;
		pathSort(array, *it, row_top, nPaths, NULL);
	}

	// run initial checker
	score = 0;
	settingToResample = NULL;
	pathLAChecker(array, path, path, 0, k, score, settingToResample, NULL);
	minCountCheck(array, c, score, settingToResample, NULL);
	cout << "Score: " << score << endl;

	for (int iter = 0; iter < iters && score > 0; iter++) {

		struct timespec start;
		struct timespec finish;
		float elapsedTime;

		// ensure we recieved an actual setting
		if (settingToResample == NULL) {
			cout << "No resampleable setting was found" << endl;
			break;
		}

		// get factors to resample (all factors in constraint group if one exists)
		conGroup = groupingInfo[settingToResample->factor_i]->conGroup;
		if (conGroup == NULL) {
			// get factor to resample
			factor_i = settingToResample->factor_i;

			// resample locating array
			for (int row_i = row_top; row_i < rows; row_i++) {
				oldLevels[row_i][factor_i] = levelMatrix[row_i][factor_i];
				if (rand() % 100 < 100) {
					levelMatrix[row_i][factor_i] = rand() % groupingInfo[factor_i]->levels;
				}
			}

			// repopulate columns of CS matrix
			for (int level_i = 0; level_i < groupingInfo[factor_i]->levels; level_i++) {
				repopulateColumns(factor_i, level_i, row_top, rows - row_top);
			}
		} else if (nCustomFactors > 0) {
			for (int row_i = row_top; row_i < rows; row_i++) {
				for (int factor_i = 0; factor_i < nCustomFactors; factor_i++) {
					oldLevels[row_i][customFactorIndeces[factor_i]] = levelMatrix[row_i][customFactorIndeces[factor_i]];
				}
				if (rand() % 100 < 100) {
					for (int factor_i = 0; factor_i < nCustomFactors; factor_i++) {
						levelMatrix[row_i][customFactorIndeces[factor_i]] = rand() % groupingInfo[customFactorIndeces[factor_i]]->levels;
					}
				}
			}

			// repopulate columns of CS matrix for every factor in constraint group and for each level
			for (int factor_i = 0; factor_i < nCustomFactors; factor_i++) {
				for (int level_i = 0; level_i < groupingInfo[customFactorIndeces[factor_i]]->levels; level_i++) {
					repopulateColumns(customFactorIndeces[factor_i], level_i, row_top, rows - row_top);
				}
			}
		} else {
			for (int row_i = row_top; row_i < rows; row_i++) {
				for (int factor_i = 0; factor_i < conGroup->factors; factor_i++) {
					oldLevels[row_i][conGroup->factorIndeces[factor_i]] = levelMatrix[row_i][conGroup->factorIndeces[factor_i]];
				}
				if (rand() % 100 < 100) {
					conGroup->randPopulateLevelRow(levelMatrix[row_i]);
				}
			}

			// repopulate columns of CS matrix for every factor in constraint group and for each level
			for (int factor_i = 0; factor_i < conGroup->factors; factor_i++) {
				for (int level_i = 0; level_i < groupingInfo[conGroup->factorIndeces[factor_i]]->levels; level_i++) {
					repopulateColumns(conGroup->factorIndeces[factor_i], level_i, row_top, rows - row_top);
				}
			}
		}


		// grab initial time
		current_utc_time( &start);
		// sort paths and recheck score
		for (std::list<Path*>::iterator it = pathList->begin(); it != pathList->end(); it++) {
			nPaths = 0;
			pathSort(array, *it, row_top, nPaths, NULL);
		}
		// check current time
		current_utc_time( &finish);
		// get elapsed seconds
		elapsedTime = (finish.tv_sec - start.tv_sec);
		// add elapsed nanoseconds
		elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
//		cout << "Elapsed After Sort: " << elapsedTime << endl;

		newScore = 0;
		FactorSetting *newSettingToResample = NULL;

		// grab initial time
		current_utc_time( &start);
		pathLAChecker(array, path, path, 0, k, newScore, newSettingToResample, NULL);
		minCountCheck(array, c, newScore, newSettingToResample, NULL);
		// check current time
		current_utc_time( &finish);
		// get elapsed seconds
		elapsedTime = (finish.tv_sec - start.tv_sec);
		// add elapsed nanoseconds
		elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
//		cout << "Elapsed After Checker: " << elapsedTime << endl;

		if (newScore <= score) { // add "|| true" to cause every change to be implemented, not just improving changes
			cout << "Rows: " << rows << " Iter: " << iter << ": ";
			cout << newScore << ": \t" << score << " \tAccepted " << endl;
			settingToResample = newSettingToResample;
			score = newScore;
		} else {
			cout << "Rows: " << rows << " Iter: " << iter << ": ";
			cout << score << " \tMaintained " << endl;

			if (conGroup == NULL) {
				// get factor to resample
				factor_i = settingToResample->factor_i;

				// rollback the change
				for (int row_i = row_top; row_i < rows; row_i++) {
					levelMatrix[row_i][factor_i] = oldLevels[row_i][factor_i];
				}

				// repopulate columns of CS matrix
				for (int level_i = 0; level_i < groupingInfo[factor_i]->levels; level_i++) {
					repopulateColumns(factor_i, level_i, row_top, rows - row_top);
				}
			} else if (nCustomFactors > 0) {
				for (int row_i = row_top; row_i < rows; row_i++) {
					for (int factor_i = 0; factor_i < nCustomFactors; factor_i++) {
						levelMatrix[row_i][customFactorIndeces[factor_i]] = oldLevels[row_i][customFactorIndeces[factor_i]];
					}
				}

				// repopulate columns of CS matrix for every factor in constraint group and for each level
				for (int factor_i = 0; factor_i < nCustomFactors; factor_i++) {
					for (int level_i = 0; level_i < groupingInfo[customFactorIndeces[factor_i]]->levels; level_i++) {
						repopulateColumns(customFactorIndeces[factor_i], level_i, row_top, rows - row_top);
					}
				}
			} else {
				// rollback the change
				for (int row_i = row_top; row_i < rows; row_i++) {
					for (int factor_i = 0; factor_i < conGroup->factors; factor_i++) {
						levelMatrix[row_i][conGroup->factorIndeces[factor_i]] = oldLevels[row_i][conGroup->factorIndeces[factor_i]];
					}
				}

				// repopulate columns of CS matrix for every factor in constraint group and for each level
				for (int factor_i = 0; factor_i < conGroup->factors; factor_i++) {
					for (int level_i = 0; level_i < groupingInfo[conGroup->factorIndeces[factor_i]]->levels; level_i++) {
						repopulateColumns(conGroup->factorIndeces[factor_i], level_i, row_top, rows - row_top);
					}
				}
			}
//			cout << score << ": \t" << newScore << " \tRejected" << endl;
		}
	}

	// deallocate all memory
	for (int row_i = 0; row_i < rows; row_i++) {
		delete[] oldLevels[row_i];
	}
	delete[] oldLevels;
	delete[] customFactorIndeces;

	// perform one last sort and save final unfinished paths to list
	list <Path*>newPathList(*pathList);
	pathList->clear();

	// sort paths and recheck score
	for (std::list<Path*>::iterator it = newPathList.begin(); it != newPathList.end(); it++) {
		nPaths = 0;
		pathSort(array, *it, row_top, nPaths, pathList);
	}
	score = 0;
	settingToResample = NULL;
	pathLAChecker(array, path, path, 0, k, score, settingToResample, NULL);
	minCountCheck(array, c, score, settingToResample, NULL);

}

// LEGACY
void CSMatrix::randomizeRows(CSCol **backupArray, CSCol **array, long long int &csScore, int row_top, int row_len) {
	int cols = getCols();
	long long int newCsScore;
	int factor_i, resampleFactor;

	FactorSetting *settingToResample = NULL;

	char **levelMatrix = locatingArray->getLevelMatrix();

	char *oldLevels = new char[rows];

	smartSort(array, row_top);
	csScore = getArrayScore(array);

	cout << "Score: " << csScore << endl;

	for (int iter = 0; iter < 1000; ) {

		if (csScore <= 0) break;

		for (int col_i = 0; col_i < cols - 1; col_i++) {

			// check if the streak ended
			if (compare(array[col_i], array[col_i + 1], 0, rows) == 0) {

				// copy to backup array
				memcpy(backupArray, array, sizeof(CSCol*) * cols); // backup array

				while (iter < 1000) {
					iter++;

					resampleFactor = rand() % (array[col_i]->factors + array[col_i + 1]->factors);

					if (resampleFactor < array[col_i]->factors) {
						factor_i = array[col_i]->setting[resampleFactor].factor_i;
					} else {
						factor_i = array[col_i + 1]->setting[resampleFactor - array[col_i]->factors].factor_i;
					}

					for (int row_i = row_top; row_i < row_top + row_len; row_i++) {
						oldLevels[row_i] = levelMatrix[row_i][factor_i];
						if (rand() % 100 < 100) {
							levelMatrix[row_i][factor_i] = rand() % groupingInfo[factor_i]->levels;
						}
					}

					for (int level_i = 0; level_i < groupingInfo[factor_i]->levels; level_i++) {
						repopulateColumns(factor_i, level_i, row_top, row_len);
					}

					smartSort(array, row_top);
					newCsScore = getArrayScore(array);

					float likelihood = 10 / pow((double)newCsScore / (double)csScore, 10);
					cout << "Rows: " << rows << " Iter: " << iter << ": ";
					if (newCsScore <= csScore) {// || rand() % 100 < likelihood) {
						cout << newCsScore << ": \t" << csScore << " \tAccepted " << likelihood << "%" << endl;
						csScore = newCsScore;
						break;
					} else {
						// rollback the change
						for (int row_i = row_top; row_i < row_top + row_len; row_i++) {
							levelMatrix[row_i][factor_i] = oldLevels[row_i];
						}

						for (int level_i = 0; level_i < groupingInfo[factor_i]->levels; level_i++) {
							repopulateColumns(factor_i, level_i, row_top, row_len);
						}

						// restore order from backup array
						memcpy(array, backupArray, sizeof(CSCol*) * cols); // backup array

						cout << csScore << ": \t" << newCsScore << " \tRejected" << endl;
					}
				}

				break;

			} else if (compare(array[col_i], array[col_i + 1], 0, rows) > 0) {
				cout << "Mistake in array" << endl;
			}
		}

	}

	delete[] oldLevels;
}

void CSMatrix::repopulateColumns(int setFactor_i, int setLevel_i, int row_top, int row_len) {

	int lastCol_i = -1;
	repopulateColumns(setFactor_i, setLevel_i, locatingArray->getFactors() - 1, locatingArray->getT(),
		mapping, locatingArray->getLevelMatrix(), lastCol_i, row_top, row_len);
}

void CSMatrix::repopulateColumns(int setFactor_i, int setLevel_i, int maxFactor_i, int t,
		Mapping *mapping, char **levelMatrix, int &lastCol_i, int row_top, int row_len) {

	if (setFactor_i > maxFactor_i && mapping->mappedTo != lastCol_i) {
		populateColumnData(data->at(mapping->mappedTo), levelMatrix, row_top, row_len);

		lastCol_i = mapping->mappedTo;
	}

	if (t == 0) return;

	int minFactor_i;
	int minLevel_i, maxLevel_i;

	if (setFactor_i > maxFactor_i) {
		minFactor_i = t - 1;
	} else {
		minFactor_i = setFactor_i;

		// if this is the last interaction, make sure setFactor shows up
		if (t == 1) maxFactor_i = minFactor_i;
	}

	for (int factor_i = minFactor_i; factor_i <= maxFactor_i; factor_i++) {

		if (factor_i == setFactor_i) {
			repopulateColumns(setFactor_i, setLevel_i, factor_i - 1, t - 1,
				mapping->mapping[factorLevelMap[factor_i][setLevel_i]], levelMatrix, lastCol_i, row_top, row_len);
		} else {
			for (int level_i = 0; level_i < groupingInfo[factor_i]->levels; level_i++) {
				repopulateColumns(setFactor_i, setLevel_i, factor_i - 1, t - 1,
					mapping->mapping[factorLevelMap[factor_i][level_i]], levelMatrix, lastCol_i, row_top, row_len);
			}
		}

	}
}

int CSMatrix::getRows() {
	return rows;
}

int CSMatrix::getCols() {
	return data->size();
}

CSCol *CSMatrix::getCol(int col_i) {
	return data->at(col_i);
}

float CSMatrix::getDistanceToCol(int col_i, float *residuals) {
	float distanceSum = 0, subtractionResult = 0;
	CSCol *csCol = data->at(col_i);

	for (int row_i = 0; row_i < rows; row_i++) {
		subtractionResult = csCol->dataP[row_i] - residuals[row_i];
		distanceSum += subtractionResult * subtractionResult;
	}

	return distanceSum;
}

float CSMatrix::getProductWithCol(int col_i, float *residuals) {
	float dotSum = 0;
	CSCol *csCol = data->at(col_i);

	for (int row_i = 0; row_i < rows; row_i++) {
		dotSum += csCol->dataP[row_i] * residuals[row_i];
	}

	return abs(dotSum);
}

string CSMatrix::getColName(CSCol *csCol) {

	ostringstream colName;

	if (csCol->factors == 0) colName << "INTERCEPT";

	// add all factor strings
	for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
		colName << (setting_i == 0 ? "" : " & ") << getFactorString(csCol->setting[setting_i]);
	}

	return colName.str();

}

string CSMatrix::getFactorString(FactorSetting setting) {
	ostringstream factorString;

	factorString << factorData->getFactorName(setting.factor_i) << "=";

	if (setting.grouped) {
		factorString << "GROUP(";
	}

	for (int index = setting.index; index < setting.index + setting.levelsInGroup; index++) {
		if (index != setting.index) factorString << "|";
		factorString << getFactorLevelName(setting.factor_i, index);
	}

	if (setting.grouped) {
		factorString << ")";
	}

	return factorString.str();
}

string CSMatrix::getFactorLevelName(int factor_i, int level_i) {
	return factorData->getFactorLevelName(factor_i, level_i);
}

void CSMatrix::addOneWayInteraction(int factor_i, char level_i, char **levelMatrix, vector <float>&sumOfSquares) {

	// create new column for CS Matrix
	CSCol *csCol = new CSCol;
	int sum = 0;

	// assign the headers
	csCol->factors = 1; // 1 contributing factor (1-way interaction)
	csCol->setting = new FactorSetting[1];
	csCol->setting[0].grouped = false;
	csCol->setting[0].factor_i = factor_i;	// single factor
	csCol->setting[0].index = level_i;		// level of single factor
	csCol->setting[0].levelsInGroup = 1;

	// populate every row
	for (int row_i = 0; row_i < rows; row_i++) {
		addRow(csCol);

		// 1 or -1 depending on if the factor levels matched
		if (level_i == levelMatrix[row_i][factor_i])
			csCol->dataP[row_i] = ENTRY_A;
		else
			csCol->dataP[row_i] = ENTRY_B;

		// add to sum of squares
		sum += csCol->dataP[row_i] * csCol->dataP[row_i];

	}

	// push into vector
	data->push_back(csCol);
	sumOfSquares.push_back(sum);

}

void CSMatrix::print() {
	cout << "CSMatrix:" << endl;
	for (int col_i = 0; col_i < getCols(); col_i++) {
		cout << getColName(data->at(col_i)) << "\t";
	}
	cout << endl;

	for (int row_i = 0; row_i < rows; row_i++) {

		for (int col_i = 0; col_i < getCols(); col_i++) {
			cout << data->at(col_i)->dataP[row_i] << "\t";
		}

		cout << endl;
	}

	// verify that the column mappings are correct
	for (int col_i = 0; col_i < getCols(); col_i++) {
		if (col_i != getColIndex(data->at(col_i))) {
			cout << "Invalid Mapping!!" << endl;
		}
	}

}

void CSMatrix::countOccurrences(CSCol *csCol, Occurrence *occurrence, int minSetting_i, float magnitude, float rSquared) {
         if (occurrence->factorList_n == csCol->factors) {
		occurrence->rSquaredContribution += rSquared;
	}
	if (occurrence->list == NULL) return;

	for (int setting_i = minSetting_i; setting_i < csCol->factors; setting_i++) {
		occurrence->list[csCol->setting[setting_i].factor_i].count++;
		occurrence->list[csCol->setting[setting_i].factor_i].magnitude += abs(magnitude);

		countOccurrences(csCol, &occurrence->list[csCol->setting[setting_i].factor_i], setting_i + 1, magnitude, rSquared);
	}
}

int CSMatrix::getColIndex(CSCol *csCol) {

	Mapping *m = mapping;
	for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
		if (m == NULL) return -1;
		m = m->mapping[factorLevelMap[csCol->setting[setting_i].factor_i][csCol->setting[setting_i].index]];
	}
	return m->mappedTo;

}

void CSMatrix::swapColumns(CSCol **array, int col_i1, int col_i2) {
	CSCol *tempCol = array[col_i1];
	array[col_i1] = array[col_i2];
	array[col_i2] = tempCol;
}

void CSMatrix::swapRows(CSCol **array, int row_i1, int row_i2) {
	float tempData;
	for (int col_i = 0; col_i < getCols(); col_i++) {
		tempData = array[col_i]->dataP[row_i1];
		array[col_i]->dataP[row_i1] = array[col_i]->dataP[row_i2];
		array[col_i]->dataP[row_i2] = tempData;
	}

	char **levelMatrix = locatingArray->getLevelMatrix();
	char *tempLevel = levelMatrix[row_i1];
	levelMatrix[row_i1] = levelMatrix[row_i2];
	levelMatrix[row_i2] = tempLevel;
}

// sort the array, given that some rows are already sorted
void CSMatrix::smartSort(CSCol **array, int sortedRows) {

	int streakMin, streakMax;
	streakMin = 0;

	for (int col_i = 1; col_i < getCols(); col_i++) {
		// check if the streak ended
		if (compare(array[col_i - 1], array[col_i], 0, sortedRows) < 0) {
			streakMax = col_i - 1;
			if (streakMin < streakMax) {
				rowSort(array, streakMin, streakMax, sortedRows, rows - sortedRows);
			}
			streakMin = col_i;
		}
	}

	// add the final streak
	streakMax = getCols() - 1;
	if (streakMin < streakMax) {
		rowSort(array, streakMin, streakMax, sortedRows, rows - sortedRows);
	}

}

void CSMatrix::rowSort(CSCol **array, int min, int max, int row_i, int row_len) {

	if (min >= max || row_len <= 0) return;

	int tempMin = min - 1;
	int tempMax = max + 1;

	while (true) {
		while (tempMin < max && array[tempMin + 1]->dataP[row_i] == ENTRY_A) tempMin++;
		while (tempMax > min && array[tempMax - 1]->dataP[row_i] == ENTRY_B) tempMax--;

		if (tempMax - 1 > tempMin + 1) {
			swapColumns(array, tempMin + 1, tempMax - 1);
		} else {
			break;
		}
	}

	rowSort(array, min, tempMin, row_i + 1, row_len - 1);
	rowSort(array, tempMax, max, row_i + 1, row_len - 1);

}

int CSMatrix::sortByCoverable(CSCol **array, int min, int max) {

	if (min >= max) return min;

	int tempMin = min - 1;
	int tempMax = max + 1;

	while (true) {
		while (tempMin < max && !array[tempMin + 1]->coverable) tempMin++;
		while (tempMax > min && array[tempMax - 1]->coverable) tempMax--;

		if (tempMax - 1 > tempMin + 1) {
			swapColumns(array, tempMin + 1, tempMax - 1);
		} else {
			break;
		}
	}

	// verification
	for (int col_i = min; col_i <= tempMin; col_i++) {
		if (array[col_i]->coverable) cout << "mistake" << endl;
	}

	return tempMin + 1;

}

int CSMatrix::sortByTWayInteraction(CSCol **array, int min, int max) {

	if (min >= max) return min;

	int t = locatingArray->getT();
	int tempMin = min - 1;
	int tempMax = max + 1;

	while (true) {
		while (tempMin < max && array[tempMin + 1]->factors < t) tempMin++;
		while (tempMax > min && array[tempMax - 1]->factors == t) tempMax--;

		if (tempMax - 1 > tempMin + 1) {
			swapColumns(array, tempMin + 1, tempMax - 1);
		} else {
			break;
		}
	}

	// verification
	for (int col_i = min; col_i <= tempMin; col_i++) {
		if (array[col_i]->factors == t) cout << "mistake" << endl;
	}

	return tempMin + 1;

}

void CSMatrix::pathSort(CSCol **array, Path *path, int row_i, int &nPaths, list <Path*>*pathList) {
	if (path->min == path->max) {
		deletePath(path->entryA);
		deletePath(path->entryB);
		path->entryA = NULL;
		path->entryB = NULL;
		return;
	} else if (row_i >= rows) {
		// add to list
		if (pathList != NULL) pathList->push_front(path);

		return;
	}

	int tempMin = path->min - 1;
	int tempMax = path->max + 1;

	while (true) {
		while (tempMin < path->max && array[tempMin + 1]->dataP[row_i] == ENTRY_A) tempMin++;
		while (tempMax > path->min && array[tempMax - 1]->dataP[row_i] == ENTRY_B) tempMax--;

		if (tempMax - 1 > tempMin + 1) {
			swapColumns(array, tempMin + 1, tempMax - 1);
		} else {
			break;
		}
	}

	// verification
	for (int col_i = path->min; col_i <= tempMin; col_i++) {
		if (array[col_i]->dataP[row_i] != ENTRY_A) cout << "mistake" << endl;
	}
	for (int col_i = tempMax; col_i <= path->max; col_i++) {
		if (array[col_i]->dataP[row_i] != ENTRY_B) cout << "mistake" << endl;
	}
	if (tempMin != tempMax - 1) cout << "mistake" << endl;

	if (path->min <= tempMin) {
		nPaths++;
		// allocate memory if none exists
		if (path->entryA == NULL) {
			path->entryA = new Path;
			path->entryA->entryA = NULL;
			path->entryA->entryB = NULL;
		}

		// populate path for entryA
		path->entryA->min = path->min;
		path->entryA->max = tempMin;

		// sort path for entryA
		pathSort(array, path->entryA, row_i + 1, nPaths, pathList);
	} else {
		// delete unnecessary path for entryA
		deletePath(path->entryA);
		path->entryA = NULL;
	}

	if (tempMax <= path->max) {
		nPaths++;
		// allocate memory if none exists
		if (path->entryB == NULL) {
			path->entryB = new Path;
			path->entryB->entryA = NULL;
			path->entryB->entryB = NULL;
		}

		// populate path for entryB
		path->entryB->min = tempMax;
		path->entryB->max = path->max;

		// sort path for entryB
		pathSort(array, path->entryB, row_i + 1, nPaths, pathList);
	} else {
		deletePath(path->entryB);
		path->entryB = NULL;
	}
}

void CSMatrix::deletePath(Path *path) {
	if (path != NULL) {
		deletePath(path->entryA);
		deletePath(path->entryB);
		delete path;
	}
}

// Locating Array Checker
void CSMatrix::pathLAChecker(CSCol **array, Path *pathA, Path *pathB, int row_i, int k,
		long long int &score, FactorSetting *&settingToResample, long long int *rowContributions) {
	if (k == 0 || pathA == NULL || pathB == NULL || pathA->min == pathB->max) {
		return;
	} else if (row_i == rows) {
//		cout << "Issue " << getColName(array[pathA->min]) << " vs " << getColName(array[pathB->max]) << endl;
		if (pathA == pathB) {
			score += (long long int)k * ((long long int)(pathA->max - pathA->min + 1) * (long long int)(pathA->max - pathA->min)) / 2;
		} else {
			score += (long long int)k * (long long int)(pathA->max - pathA->min + 1) * (long long int)(pathB->max - pathB->min + 1);
		}

		// set a setting to resample
		if (settingToResample == NULL) {
			int columnToResample = -1;

			/*
			int offset;
			do {
				offset = rand() % (pathA->max - pathA->min + 1 + pathB->max - pathB->min + 1);
				if (offset <= pathA->max - pathA->min) {
					columnToResample = pathA->min + offset;
				} else {
					offset -= (pathA->max - pathA->min + 1);
					columnToResample = pathB->min + offset;
				}
			} while (array[columnToResample]->factors <= 0);
			*/

			// find a pair of columns that is distinguishable
			for (int i_a = pathA->min;
					i_a <= pathA->max && columnToResample == -1; i_a++) {
				for (int i_b = (pathA == pathB ? i_a + 1 : 0);
						i_b <= pathB->max && columnToResample == -1; i_b++) {
					if (checkDistinguishable(array[i_a], array[i_b])) {
						if (rand() % 2 && array[i_a]->factors > 0) {
							columnToResample = i_a;
						} else if (array[i_b]->factors > 0) {
							columnToResample = i_b;
						}
					} else {
//						score--;
//						cout << "Not distinguishable: " << getColName(array[i_a]) << " vs " << getColName(array[i_b]) << endl;
					}
				}
			}

			// randomly choose a setting in the column to resample
			if (columnToResample != -1) {
				settingToResample = &array[columnToResample]->setting[rand() % array[columnToResample]->factors];
			}
		}

		return;
	}

	Path *pathAentryA, *pathAentryB, *pathBentryA, *pathBentryB;

	if (pathA->min == pathA->max) {
		pathAentryA = (array[pathA->min]->dataP[row_i] == ENTRY_A ? pathA : NULL);
		pathAentryB = (array[pathA->min]->dataP[row_i] == ENTRY_B ? pathA : NULL);
	} else {
		pathAentryA = pathA->entryA;
		pathAentryB = pathA->entryB;
	}

	if (pathB->min == pathB->max) {
		pathBentryA = (array[pathB->min]->dataP[row_i] == ENTRY_A ? pathB : NULL);
		pathBentryB = (array[pathB->min]->dataP[row_i] == ENTRY_B ? pathB : NULL);
	} else {
		pathBentryA = pathB->entryA;
		pathBentryB = pathB->entryB;
	}

	pathLAChecker(array, pathAentryA, pathBentryA, row_i + 1, k, score, settingToResample, rowContributions);
	pathLAChecker(array, pathAentryB, pathBentryB, row_i + 1, k, score, settingToResample, rowContributions);
	pathLAChecker(array, pathAentryA, pathBentryB, row_i + 1, k - 1, score, settingToResample, rowContributions);

	// add row contributions
	if (rowContributions != NULL && pathAentryA != NULL && pathBentryB != NULL) {
		rowContributions[row_i] += (pathAentryA->max - pathAentryA->min + 1) * (pathBentryB->max - pathBentryB->min + 1);
	}

	if (pathA != pathB) {
		pathLAChecker(array, pathAentryB, pathBentryA, row_i + 1, k - 1, score, settingToResample, rowContributions);

		// add row contributions
		if (rowContributions != NULL && pathAentryB != NULL && pathBentryA != NULL) {
			rowContributions[row_i] += (pathAentryB->max - pathAentryB->min + 1) * (pathBentryA->max - pathBentryA->min + 1);
		}
	}
}

// Detecting Array Checker
void CSMatrix::pathDAChecker(CSCol **array, Path *pathA, Path *pathB, int row_i, int k,
		long long int &score, FactorSetting *&settingToResample, long long int *rowContributions) {
	if (k == 0 || pathA == NULL || pathB == NULL || pathA->min == pathB->max) {
		return;
	} else if (row_i == rows) {
//		cout << getColName(array[pathA->min]) << " vs " << getColName(array[pathB->max]) << endl;
		if (pathA == pathB) {
			score += (long long int)k * ((long long int)(pathA->max - pathA->min + 1) * (long long int)(pathA->max - pathA->min)) / 2;
		} else {
			score += (long long int)k * (long long int)(pathA->max - pathA->min + 1) * (long long int)(pathB->max - pathB->min + 1);
		}

		// set a setting to resample
		if (settingToResample == NULL) {
			int columnToResample;
			int offset;
			do {
				offset = rand() % (pathA->max - pathA->min + 1 + pathB->max - pathB->min + 1);
				if (offset <= pathA->max - pathA->min) {
					columnToResample = pathA->min + offset;
				} else {
					offset -= (pathA->max - pathA->min + 1);
					columnToResample = pathB->min + offset;
				}
			} while (array[columnToResample]->factors <= 0);

			settingToResample = &array[columnToResample]->setting[rand() % array[columnToResample]->factors];
		}

		return;
	}

	Path *pathAentryA, *pathAentryB, *pathBentryA, *pathBentryB;

	if (pathA->min == pathA->max) {
		pathAentryA = (array[pathA->min]->dataP[row_i] == ENTRY_A ? pathA : NULL);
		pathAentryB = (array[pathA->min]->dataP[row_i] == ENTRY_B ? pathA : NULL);
	} else {
		pathAentryA = pathA->entryA;
		pathAentryB = pathA->entryB;
	}

	if (pathB->min == pathB->max) {
		pathBentryA = (array[pathB->min]->dataP[row_i] == ENTRY_A ? pathB : NULL);
		pathBentryB = (array[pathB->min]->dataP[row_i] == ENTRY_B ? pathB : NULL);
	} else {
		pathBentryA = pathB->entryA;
		pathBentryB = pathB->entryB;
	}

	pathDAChecker(array, pathAentryA, pathBentryA, row_i + 1, k, score, settingToResample, rowContributions);
	pathDAChecker(array, pathAentryB, pathBentryB, row_i + 1, k, score, settingToResample, rowContributions);
	pathDAChecker(array, pathAentryB, pathBentryA, row_i + 1, k, score, settingToResample, rowContributions);
	pathDAChecker(array, pathAentryA, pathBentryB, row_i + 1, k - 1, score, settingToResample, rowContributions);

	// add row contributions
	if (rowContributions != NULL && pathAentryA != NULL && pathBentryB != NULL) {
		rowContributions[row_i] += (pathAentryA->max - pathAentryA->min + 1) * (pathBentryB->max - pathBentryB->min + 1);
	}

}

void CSMatrix::quickSort(CSCol **array, int min, int max, int row_top, int row_len) {
	if (min == max) return;

	int tempMin = min;
	int tempMax = max;
	CSCol *midCol = array[(min + max) / 2];

	while (true) {

		while (tempMin <= tempMax && tempMin < max && compare(array[tempMin], midCol, row_top, row_len) < 0) {
			tempMin++;
		}

		while (tempMin <= tempMax && tempMax > min && compare(array[tempMax], midCol, row_top, row_len) > 0) {
			tempMax--;
		}

		if (tempMin == tempMax) {
			if (tempMin == min) tempMin++;
			else tempMax--;
			break;
		} else if (tempMin <= tempMax) {
			swapColumns(array, tempMin, tempMax);
			tempMin++;
			tempMax--;
		} else {
			break;
		}

	}

	quickSort(array, min, tempMax, row_top, row_len);
	quickSort(array, tempMin, max, row_top, row_len);
}

int CSMatrix::compare(CSCol *csCol1, CSCol *csCol2, int row_top, int row_len) {
	int length = row_len * sizeof(float);
	return memcmp(&csCol1->dataP[row_top], &csCol2->dataP[row_top], length);
}

// LEGACY (to be removed later)
//cout << "Checking advanced for " << k << " differences" << endl;
//score = checkAdvanced(array, k, 0, getCols() - 1, 0, rows, settingToResample);
//cout << "Advanced Elapsed: " << elapsedTime << endl;
long CSMatrix::checkAdvanced(CSCol **array, int k, int min, int max, int row_top, int row_len, FactorSetting *&settingToResample) {

	// no more differences to find or only 1 column
	if (k <= 0 || min >= max) {
		return 0;
	} else if (k > row_len) {

//		cout << "Issue " << getColName(array[min]) << " vs " << getColName(array[max]) << endl;
//		cout << "Looking for " << k << " differences, but there are " << row_len << " rows left " << min << " vs " << max << endl;

		CSCol *csCol1 = array[min];
		CSCol *csCol2 = array[max];

		if (settingToResample == NULL) {
			// choose random factor+level for resampling
			CSCol *csColToResample = NULL;

			if (csCol1->factors != 0 && csCol2->factors != 0) {
				// neither has 0 factors, choose randomly
				if (rand() % 2) {
					csColToResample = csCol1;
				} else {
					csColToResample = csCol2;
				}
			} else if (csCol1->factors != 0) {
				// csCol1 has factors so choose it (csCol2 has no factors)
				csColToResample = csCol1;
			} else if (csCol2->factors != 0) {
				// csCol2 has factors so choose it (csCol1 has no factors)
				csColToResample = csCol2;
			}

			// chose a random factor setting to resample
			if (csColToResample != NULL) {
				settingToResample = &csColToResample->setting[rand() % csColToResample->factors];
			}
		}

		int cols = max - min + 1;
		return (long) pow((double) (cols * cols), (double) k);
	} else {

		long score = 0;

		rowSort(array, min, max, row_top, 1);

		int divide;
		for (divide = min; divide < max; divide++) {
			if (compare(array[divide], array[divide + 1], row_top, 1) != 0) {
				break;
			}
		}
		for (int col_i = divide + 1; col_i < max; col_i++) {
			if (compare(array[col_i], array[col_i + 1], row_top, 1) != 0) {
				cout << "Integrity issue" << endl;
			}
		}

		// check the left side
		score += checkAdvanced(array, k, min, divide, row_top + 1, row_len - 1, settingToResample);
		// check the right side
		score += checkAdvanced(array, k, divide + 1, max, row_top + 1, row_len - 1, settingToResample);
		// check both
		score += checkAdvanced(array, k - 1, min, max, row_top + 1, row_len - 1, settingToResample);

		return score;
	}
}

// FUNCTIONS FOR 'fixla'

void CSMatrix::addRow(CSCol **array, char *levelRow) {
	rows++;

	locatingArray->addLevelRow(levelRow);

	// add a row to each column of the CS matrix and populate
	for (int col_i = 0; col_i < getCols(); col_i++) {
		addRow(array[col_i]);
		populateColumnData(array[col_i], locatingArray->getLevelMatrix(), rows - 1, 1);
	}
}

void CSMatrix::remRow(CSCol **array) {
	rows--;

	char *levelRow = locatingArray->remLevelRow();
	delete levelRow;

	// remove a row from each column of the CS matrix
	for (int col_i = 0; col_i < getCols(); col_i++) {
		remRow(array[col_i]);
	}
}

void CSMatrix::resizeArray(CSCol **array, int newRows) {
	int factors = locatingArray->getFactors();
	int nConGroups = locatingArray->getNConGroups();
	ConstraintGroup **conGroups = locatingArray->getConGroups();

	// increase size as needed
	while (newRows > rows) {
		// allocate memory for new row of locating array
		char *levelRow = new char[factors];

		// generate random row
		for (int factor_i = 0; factor_i < factors; factor_i++) {
			levelRow[factor_i] = rand() % groupingInfo[factor_i]->levels;
		}

		// resample constraint groups
		for (int conGroup_i = 0; conGroup_i < nConGroups; conGroup_i++) {
			conGroups[conGroup_i]->randPopulateLevelRow(levelRow);
		}

		addRow(array, levelRow);
	}

	// decrease size as needed
	while (newRows < rows && rows > 1) {
		remRow(array);
	}
}

bool CSMatrix::checkColumnCoverability(CSCol *csCol) {

	int factors = locatingArray->getFactors();
	int nConGroups = locatingArray->getNConGroups();
	ConstraintGroup **conGroups = locatingArray->getConGroups();

	char *requireLevelRow = new char[factors];
	char *avoidLevelRow = new char[factors];
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		requireLevelRow[factor_i] = -1;
		avoidLevelRow[factor_i] = -1;
	}

	// set to require csCol1 settings
	for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
		requireLevelRow[csCol->setting[setting_i].factor_i] = csCol->setting[setting_i].index;
	}

	// check if possible for every constraint group
	bool satisfiableInGroups = true;
	for (int conGroup_i = 0; conGroup_i < nConGroups; conGroup_i++) {
		if (!conGroups[conGroup_i]->satisfiableInGroupLA(requireLevelRow, avoidLevelRow)) {
			satisfiableInGroups = false;
			break;
		}
	}

	/* for factors not in a constraint group, simply set their settings in csCol1 to be ENTRY_A,
	 and the check passes because the interaction cannot be only */

	// deallocate memory
	delete[] requireLevelRow;
	delete[] avoidLevelRow;

	return satisfiableInGroups;
}

bool CSMatrix::checkDistinguishable(CSCol *csCol1, CSCol *csCol2) {
	return checkOneWayDistinguishable(csCol1, csCol2) || checkOneWayDistinguishable(csCol2, csCol1);
}

bool CSMatrix::checkOneWayDistinguishable(CSCol *csCol1, CSCol *csCol2) {

	int factors = locatingArray->getFactors();
	int nConGroups = locatingArray->getNConGroups();
	ConstraintGroup **conGroups = locatingArray->getConGroups();

	char *requireLevelRow = new char[factors];
	char *avoidLevelRow = new char[factors];
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		requireLevelRow[factor_i] = -1;
		avoidLevelRow[factor_i] = -1;
	}

	// set to require csCol1 settings
	for (int setting_i = 0; setting_i < csCol1->factors; setting_i++) {
		requireLevelRow[csCol1->setting[setting_i].factor_i] = csCol1->setting[setting_i].index;
	}

	// at least one satisfying setting in csCol2 must be avoidable for every constraint group
	bool settingExists = false;

	// try to avoid at least one csCol2 setting
	for (int setting_i = 0; setting_i < csCol2->factors; setting_i++) {
		// we are trying to avoid "csCol2->setting[setting_i].factor_i != csCol2->setting[setting_i].index" in this iteration
		// we cannot require the setting we are trying to avoid
		if (requireLevelRow[csCol2->setting[setting_i].factor_i] != csCol2->setting[setting_i].index) {
			avoidLevelRow[csCol2->setting[setting_i].factor_i] = csCol2->setting[setting_i].index;

			// check if possible for every constraint group
			bool satisfiableInGroups = true;
			for (int conGroup_i = 0; conGroup_i < nConGroups; conGroup_i++) {
				if (!conGroups[conGroup_i]->satisfiableInGroupLA(requireLevelRow, avoidLevelRow)) {
					satisfiableInGroups = false;
					break;
				}
			}

			if (satisfiableInGroups) {
				settingExists = true;
				break;
			}

			avoidLevelRow[csCol2->setting[setting_i].factor_i] = -1;
		}
	}

	/* for factors not in a constraint group, simply set their settings in csCol1 to be ENTRY_A,
	 and the check passes because the interaction cannot be only */

	// deallocate memory
	delete[] requireLevelRow;
	delete[] avoidLevelRow;

	return settingExists;
}

void CSMatrix::randomizeArray(CSCol **array) {
	int factors = locatingArray->getFactors();
	int cols = getCols();

	// resample entire locating array
	char **levelMatrix = locatingArray->getLevelMatrix();
	for (int row_i = 0; row_i < rows; row_i++) {
		for (int factor_i = 0; factor_i < factors; factor_i++) {
			levelMatrix[row_i][factor_i] = rand() % groupingInfo[factor_i]->levels;
		}
	}

	// populate each column of the CS matrix
	for (int col_i = 0; col_i < cols; col_i++) {
		populateColumnData(array[col_i], locatingArray->getLevelMatrix(), 0, rows);
	}
}

void CSMatrix::addRowFix(CSCol **array, long long int &csScore) {
	// the maximum seconds without finding a better score
	int secsMax = 2;

	// total factors in locating array
	int factors = locatingArray->getFactors();

	long long int bestScore;
	long long int newScore;
	CSCol *csCol, *bestCol;

	int cols = getCols();

	// for backing up the order of array
	CSCol **backupArray = new CSCol*[cols];

	// allocate memory for new row of locating array
	char *levelRow = new char[factors];
	char *oldLevelRow = new char[factors];
	char *newLevelRow = new char[factors];
	char *bestLevelRow = new char[factors];

	// track if factors of new row are finalized
	bool *finalized = new bool[factors];

	// generate random row non-finalized
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		finalized[factor_i] = false;

		levelRow[factor_i] = rand() % groupingInfo[factor_i]->levels;
	}

	// add the row to locating array
	addRow(array, levelRow);
	cout << "The matrix now has " << rows << " rows" << endl;

	// smartly sort the array and score it
	smartSort(array, rows - 1);
	csScore = getArrayScore(array);

//	cout << "Score after random non-finalized row: " << csScore << endl;

	// used to track duplicate columns in CS matrix
	bool duplicate = false;
	bool lastPairMatched = false;

	// continue adding finalizing factors while CS score decreases (improves)
	while (true) {

		// find the best column index (with the best score)
		bestScore = csScore;	// the original best is the original score
		bestCol = NULL;			// no best column yet

		// backup the sorted order of the array
		memcpy(backupArray, array, sizeof(CSCol*) * cols); // backup array

		// grab initial time
		struct timespec start;
		struct timespec finish;
		float elapsedTime;
		current_utc_time( &start);

		// go through all columns of CS matrix
		for (int col_i = 0; col_i < getCols() - 1; col_i++) {

			// check current time
			current_utc_time( &finish);
			// get elapsed seconds
			elapsedTime = (finish.tv_sec - start.tv_sec);
			// add elapsed nanoseconds
			elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
			if (elapsedTime > secsMax) break;

			// check for duplicate column
			if (compare(array[col_i], array[col_i + 1], 0, rows) >= 0) {
				duplicate = true;
				lastPairMatched = true;
			} else if (lastPairMatched) {
				duplicate = true;
				lastPairMatched = false;
			} else {
				duplicate = false;
			}

			// retrieve CSCol (and duplicate) from sorted array
			csCol = array[col_i];

			// if duplicate and the column has factors to change (column is not the INTERCEPT)
			if (duplicate && csCol->dataP[rows - 1] != 1 && csCol->factors > 0) {

				// grab duplicate column
				CSCol *csDup = array[col_i + 1];

				// the goal is to make the last row of csCol be a 1
				// all factors must be changed, check if this change conflicts with finalized factors
				bool changeAllowed = true;
				int factor_i;

				// go through all relevant factors for this column
				for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
					factor_i = csCol->setting[setting_i].factor_i;

					// backup the factor level
					oldLevelRow[factor_i] = levelRow[factor_i];

					// check if finalized yet
					if (!finalized[factor_i]) {
						// update factor to match the current column
						newLevelRow[factor_i] = csCol->setting[setting_i].index +
							(rand() % csCol->setting[setting_i].levelsInGroup);
					} else {
						// ensure this will actually make this column a 1
						changeAllowed &= (newLevelRow[factor_i] >= csCol->setting[setting_i].index &&
							newLevelRow[factor_i] < csCol->setting[setting_i].index + csCol->setting[setting_i].levelsInGroup);
					}

				}

				// check if this column change is allowed
				if (changeAllowed) {
					// CAREFUL!!!
					// col_i cannot be used after smartSort and before the rollback
					// because the array is reordered. Use csCol

					// try making the changes and see if the CS score gets smaller
					newScore = csScore;
					for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
						factor_i = csCol->setting[setting_i].factor_i;

						if (levelRow[factor_i] != newLevelRow[factor_i]) {
							// update columns
							levelRow[factor_i] = newLevelRow[factor_i];
							repopulateColumns(factor_i, oldLevelRow[factor_i], rows - 1, 1);
							repopulateColumns(factor_i, newLevelRow[factor_i], rows - 1, 1);
						}
					}

					// smart sort and get new score
					smartSort(array, rows - 1);
					newScore = getArrayScore(array);

					// check if we have a better score
					if (newScore < bestScore) {
						// store the better score as the best, and save the column
						bestScore = newScore;
						bestCol = csCol;

						for (int factor_i = 0; factor_i < factors; factor_i++) {
							bestLevelRow[factor_i] = newLevelRow[factor_i];
						}
					}

					// rollback the change
					for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
						factor_i = csCol->setting[setting_i].factor_i;

						if (oldLevelRow[factor_i] != levelRow[factor_i]) {
							// rollback columns and get new score
							levelRow[factor_i] = oldLevelRow[factor_i];
							repopulateColumns(factor_i, oldLevelRow[factor_i], rows - 1, 1);
							repopulateColumns(factor_i, newLevelRow[factor_i], rows - 1, 1);
						}
					}

					// restore old order
					memcpy(array, backupArray, sizeof(CSCol*) * getCols());

				}
			}
		}

		// check if we found a better score
		if (bestScore < csScore) {
//			cout << "Found better score! " << bestScore << endl;

			// make the change to improve the score
			int factor_i;
			for (int setting_i = 0; setting_i < bestCol->factors; setting_i++) {
				factor_i = bestCol->setting[setting_i].factor_i;

				if (bestLevelRow[factor_i] != levelRow[factor_i]) {
					// update columns and get new score
					oldLevelRow[factor_i] = levelRow[factor_i];
					levelRow[factor_i] = bestLevelRow[factor_i];
					repopulateColumns(factor_i, oldLevelRow[factor_i], rows - 1, 1);
					repopulateColumns(factor_i, bestLevelRow[factor_i], rows - 1, 1);
				}

				// finalize the changes
				finalized[factor_i] = true;
			}

			// do a final smart sort and finalize the changed factors
			smartSort(array, rows - 1);
			csScore = getArrayScore(array);

		} else {
			break;
		}

	}

	cout << "Score after finalized row: " << csScore << endl;

	delete[] oldLevelRow;
	delete[] newLevelRow;
	delete[] bestLevelRow;
	delete[] backupArray;
	delete[] finalized;

}

long long int CSMatrix::getArrayScore(CSCol **array) {
	long long int streak = 0;
	long long int squaredSum = 0;

	for (int col_i = 0; col_i < getCols() - 1; col_i++) {

		// increment streak
		streak++;

		// check if the streak ended
		if (compare(array[col_i], array[col_i + 1], 0, rows) < 0) {
			squaredSum += streak * streak;
			streak = 0;
		} else if (compare(array[col_i], array[col_i + 1], 0, rows) > 0) {
			cout << "Mistake in CS matrix at column: " << col_i << endl;
		} else {
			if (checkDistinguishable(array[col_i], array[col_i + 1])) {
//				cout << "Duplicates found: " << getColName(array[col_i]) << " vs " << getColName(array[col_i + 1]) << endl;
			}
		}
	}

	// add the final streak
	streak++;
	squaredSum += streak * streak;

	return (squaredSum - getCols());
}

long long int CSMatrix::getBruteForceArrayScore(CSCol **array, int k) {
	long long int score = 0;
	int differences;

	int indistinguishable = 0;

	for (int col_i1 = 0; col_i1 < getCols() - 1; col_i1++) {
		for (int col_i2 = col_i1 + 1; col_i2 < getCols(); col_i2++) {
			if (array[col_i1]->coverable && array[col_i2]->coverable) {
				if (!checkDistinguishable(array[col_i1], array[col_i2])) {
					cout << "Indistinguishable pair: " << getColName(array[col_i1]) << " vs " << getColName(array[col_i2]) << endl;
					indistinguishable++;
				}

				differences = 0;

				for (int row_i = 0; row_i < rows; row_i++) {
					if (array[col_i1]->dataP[row_i] != array[col_i2]->dataP[row_i]) differences++;
				}

				if (differences < k) {
					score += (k - differences);
				}
			}
		}
	}

	cout << "Indistinguishable pairs: " << indistinguishable << endl;
	cout << "Minimum score: " << (indistinguishable * k) << endl;

	return score;
}

void CSMatrix::writeResponse(string responseDir, string responseCol, int terms, float *coefficients, int *columns) {
	cout << responseDir << endl;

	// initialize all responses to 0
	float *responses = new float[rows];
	for (int row_i = 0; row_i < rows; row_i++) {
		responses[row_i] = 0;
	}

	// add contributions from all terms
	for (int term_i = 0; term_i < terms; term_i++) {

		CSCol *csCol = data->at(columns[term_i]);
		float coefficient = coefficients[term_i];

		cout << coefficient << " * " << getColName(csCol) << endl;

		for (int row_i = 0; row_i < rows; row_i++) {
			responses[row_i] += coefficient * csCol->dataP[row_i];
		}
	}

	// open a response file for writing
	string file = responseDir + "/Response.tsv";
	ofstream ofs(file.c_str());

	// write response headers
	ofs << rows << endl;
	ofs << responseCol << endl;

	// write response values
	for (int row_i = 0; row_i < rows; row_i++) {
		ofs << responses[row_i] << endl;
	}

	// close the file
	ofs.close();

}



CSMatrix::~CSMatrix() {
	CSCol *csCol;

	// delete all mappings strategically
	for (int t = locatingArray->getT(); t >= 0; t--) {
		for (int col_i = 0; col_i < getCols(); col_i++) {
			csCol = data->at(col_i);

			if (csCol->factors == t) {
				Mapping *mappingTemp = mapping;
				FactorSetting *setting;
				for (int setting_i = 0; setting_i < t; setting_i++) {
					setting = &csCol->setting[setting_i];
					mappingTemp = mappingTemp->mapping[factorLevelMap[setting->factor_i][setting->index]];
				}

				// deallocate mapping->mapping array
				if (t == 1 && groupingInfo[setting->factor_i]->grouped && setting->index > 0 &&
						groupingInfo[setting->factor_i]->levelGroups[setting->index] == groupingInfo[setting->factor_i]->levelGroups[setting->index - 1]) {
					// This column is a main effect that is grouped with other main effects, but it is not the first main effect in its group.
					// All main effects in a particular group have separate columns in the CS matrix and separate mapping structs but they share mapping->mapping arrays,
					// because the main effects (although separate) are not differentiated in any interactions that they are included in.
					// We only deallocate the mapping->mapping array for the first main effect in a group so as not to deallocate the same array more than once.
					// We therefore ignore this main effect in terms of deallocating the mapping-mapping array.
				} else if (mappingTemp->mapping != NULL) {
					delete[] mappingTemp->mapping;
				}

				// deallocate mapping struct
				delete mappingTemp;
			}
		}
	}

	// delete all csCols
	for (int col_i = 0; col_i < getCols(); col_i++) {
		csCol = data->at(col_i);

		delete[] csCol->setting;
		delete csCol;
	}
	delete data;

	// delete factor-level map
	for (int factor_i = 0; factor_i < locatingArray->getFactors(); factor_i++) {
		delete[] factorLevelMap[factor_i];
	}
	delete[] factorLevelMap;

}


// [[Rcpp::export]]
SEXP makeCSMatrix2(SEXP la){
  Rcpp::XPtr<LocatingArray> la_ptr(la);
  CSMatrix* cs = new CSMatrix(la_ptr);
  Rcpp::XPtr<CSMatrix> ptr(cs);
  return ptr;
}
