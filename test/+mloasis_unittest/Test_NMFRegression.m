classdef Test_NMFRegression < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 20-Apr-2023 21:26:42 by jjlee in repository /Users/jjlee/MATLAB-Drive/mloasis/test/+mloasis_unittest.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mloasis.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
    end
    
    methods (TestClassSetup)
        function setupNMFRegression(this)
            import mloasis.*
            this.testObj_ = NMFRegression();
        end
    end
    
    methods (TestMethodSetup)
        function setupNMFRegressionTest(this)
            this.testObj = this.testObj_;
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
