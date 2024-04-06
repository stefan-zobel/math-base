/*
 * Copyright 2018 Stefan Zobel
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
package math.trmm;

final class BlockSizes {

    static final int MR = 4; // 4
    static final int NR = 4; // 4

    static final int MC = 384; // 384
    static final int KC = 384; // 384
    static final int NC = 4096; // 4096 .. 16384

    private BlockSizes() {
        throw new AssertionError();
    }
}
