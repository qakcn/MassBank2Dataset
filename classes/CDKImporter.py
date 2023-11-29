# This script contains functions that used to import Java package CDK.
#
# Author: qakcn
# Email: qakcn@hotmail.com
# Date: 2023-11-29

# Copyright 2023 qakcn
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


if __name__ == "__main__":
    raise SystemExit("This script is not meant to be run directly")

import jpype
import jpype.imports

def start_jvm():
    if not jpype.isJVMStarted():
        jpype.startJVM("-ea", classpath = "./classes/jars/cdk-2.9.jar")

def attach_jvm():
    if not jpype.isThreadAttachedToJVM():
        jpype.attachThreadToJVM()

# pyright: reportMissingImports=false
def import_cdk():
    import org.openscience.cdk as cdk
    return cdk

def shutdown_jvm():
    jpype.shutdownJVM()


