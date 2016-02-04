/** librf main include
 * @file 
 *
 */


/**
 * \mainpage librf docs 
 * \section intro_sec Introduction
 * This is the documentation (in progress) for librf.
 *
 * The important classes are:
 * librf::RandomForest and librf::InstanceSet
 * \section example Example programs
 * \subsection rftrain  rftrain
 * Train a randomforest from a .libsvm or a .csv file.  Saves the
 * randomforest to a file.
 * \subsection rfpredict rfpredict
 * Predict probabilities using a saved random forest.
 * \subsection rfvarimport rfvarimport
 * Train a random forest and determine variable importances.
 *
 */
#ifndef _LIBRF_H_
#define _LIBRF_H_

#include "librf/random_forest.h"
#include "librf/instance_set.h"

#endif  // LIBRF_H
