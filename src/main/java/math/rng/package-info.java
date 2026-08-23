/*
 * Copyright 2013 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/**
 * Pseudo-random number generators.
 * <p>
 * Every generator implements {@link math.rng.PseudoRandom}, which documents the
 * single draws, and inherits the distribution streams of
 * {@link math.rng.PseudoRandomStream}. {@link math.rng.DefaultRng} is the entry
 * point for callers with no reason to pick a particular algorithm.
 *
 * <h2>Seeding</h2>
 *
 * Each generator offers the same three constructors:
 * <ul>
 * <li>the no-argument one seeds from {@link math.rng.SplitMix64Seed#seed()},
 * which is unpredictable and therefore not reproducible;</li>
 * <li>{@code (long seed)} and {@code (long[] seed)} are reproducible -- the
 * same seed always yields the same sequence.</li>
 * </ul>
 *
 * The array form has two readings, and which one applies depends on the
 * generator:
 * <ul>
 * <li>{@link math.rng.MersenneTwister64} consumes the whole array as seed
 * material, so a longer array carries more of it. Having nothing to seed with
 * is an error there: {@code null} raises a {@code NullPointerException} and an
 * empty array an {@code IllegalArgumentException}.</li>
 * <li>Every other generator hashes the array down to one {@code long} through
 * {@link math.rng.SplitMix64Seed#seed(long[])}, whose state is a single
 * {@code long} anyway. That mapping sends {@code null}, {@code new long[0]} and
 * {@code new long[] {0}} to the same fixed seed, so the three are
 * interchangeable and none of them is rejected.</li>
 * </ul>
 *
 * The seed a generator was built from is reported by
 * {@link math.rng.PseudoRandom#getSeed()}.
 *
 * <h2>Threads</h2>
 *
 * None of the generators is thread-safe, and sharing one is not merely a
 * quality problem: a torn multi-word state can leave a generator outside its
 * legal state space. Give each thread its own instance -- from
 * {@link math.rng.DefaultRng#newIndepPseudoRandoms(int)}, from
 * {@link math.rng.SplittablePseudoRandom#split()}, or simply by constructing
 * one per thread, which costs microseconds at most.
 * <p>
 * The streams may be run in parallel: each split draws from its own
 * generator, at the price that the values then depend on the parallelism of
 * the common {@link java.util.concurrent.ForkJoinPool} and are no longer
 * reproducible from the seed.
 */
package math.rng;
