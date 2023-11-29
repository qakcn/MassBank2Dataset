if __name__ == '__main__':
    raise SystemExit('This script is not meant to be run directly')

import jpype
import jpype.imports

# pyright: reportMissingImports=false
if not jpype.isJVMStarted():
    jpype.startJVM("-ea", classpath = './classes/jars/cdk-2.9.jar')
    
    import org.openscience.cdk as cdk
