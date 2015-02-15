//
//  MKSONLESolution.m
//  SONLE
//
//  Created by Андрей Рычков on 16.03.14.
//  Copyright (c) 2014 Andrey Rychkov. All rights reserved.
//

#import "MKSONLESolution.h"
#import "MKMatrix.h"

@implementation MKSONLESolution

- (id)initWithEpsilon:(NSNumber *)eps
        approximation:(NSArray *)approximation
    numberOfVariables:(unsigned long)number
{
    self = [super init];
    
    if (self)
    {
        self.epsilon = eps;
        self.approximation = approximation;
        self.jacobi = [[MKMatrix alloc] init];
        self.function = [NSMutableArray arrayWithCapacity:number];
        self.previousDelta = [NSMutableArray arrayWithCapacity:number];
        
        for (long i = 0; i < number; ++i)
        {
            [self.jacobi.matrix addObject:[NSMutableArray arrayWithCapacity:number]];
        }
    }
    
    return self;
}

// –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

- (void)calculateFunctionInPoint:(NSArray *)point
{
    double x1 = [point[0] doubleValue];
    double x2 = [point[1] doubleValue];
    double x3 = [point[2] doubleValue];
    double x4 = [point[3] doubleValue];
    double x5 = [point[4] doubleValue];
    double x6 = [point[5] doubleValue];
    double x7 = [point[6] doubleValue];
    double x8 = [point[7] doubleValue];
    double x9 = [point[8] doubleValue];
    double x10 = [point[9] doubleValue];
    
    self.function[0] = [NSNumber numberWithDouble:-(cos(x1 * x2) - exp(-3.0 * x3) + x4 * x5 * x5 - x6 - sinh(2.0 * x8) * x9 + 2.0 * x10 + 2.0004339741653854440)];
    self.function[1] = [NSNumber numberWithDouble:-(sin(x1 * x2) + x3 * x9 * x7 - exp(-x10 + x6) + 3.0 * x5 * x5 - x6 * (x8 + 1.0) + 10.886272036407019994)];
    self.function[2] = [NSNumber numberWithDouble:-(x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8 + x9 - x10 - 3.1361904761904761904)];
    self.function[3] = [NSNumber numberWithDouble:-(2.0 * cos(-x9 + x4) + x5 / (x3 + x1) - sin(pow(x2, 2.0)) + pow(cos(x7 * x10), 2.0) - x8 - 0.1707472705022304757)];
    self.function[4] = [NSNumber numberWithDouble:-(sin(x5) + 2.0 * x8 * (x3 + x1) - exp(-x7 * (-x10 + x6)) + 2.0 * cos(x2) - 1.0 / (x4 - x9) - 0.3685896273101277862)];
    self.function[5] = [NSNumber numberWithDouble:-(exp(x1 - x4 - x9) + x5 * x5 / x8 + cos(3.0 * x10 * x2) / 2.0 - x6 * x3 + 2.0491086016771875115)];
    self.function[6] = [NSNumber numberWithDouble:-(pow(x2, 3.0) * x7 - sin(x10 / x5 + x8) + (x1 - x6) * cos(x4) + x3 - 0.7380430076202798014)];
    self.function[7] = [NSNumber numberWithDouble:-(x5 * pow(x1 - 2.0 * x6, 2.0) - 2.0 * sin(-x9 + x3) + 1.5 * x4 - exp(x2 * x7 + x10) + 3.5668321989693809040)];
    self.function[8] = [NSNumber numberWithDouble:-(7.0 / x6 + exp(x5 + x4) - 2.0 * x2 * x8 * x10 * x7 + 3.0 * x9 - 3.0 * x1 - 8.4394734508383257499)];
    self.function[9] = [NSNumber numberWithDouble:-(x10 * x1 + x9 * x2 - x8 * x3 + sin(x4 + x5 + x6) * x7 - 0.78238095238095238096)];
}

- (void)calculateJacobiInPoint:(NSArray *)point
{
    double x1 = [point[0] doubleValue];
    double x2 = [point[1] doubleValue];
    double x3 = [point[2] doubleValue];
    double x4 = [point[3] doubleValue];
    double x5 = [point[4] doubleValue];
    double x6 = [point[5] doubleValue];
    double x7 = [point[6] doubleValue];
    double x8 = [point[7] doubleValue];
    double x9 = [point[8] doubleValue];
    double x10 = [point[9] doubleValue];
    
    self.jacobi.matrix[0][0] = [NSNumber numberWithDouble:-sin(x1 * x2) * x2];
    self.jacobi.matrix[0][1] = [NSNumber numberWithDouble:-sin(x1 * x2) * x1];
    self.jacobi.matrix[0][2] = [NSNumber numberWithDouble:3.0 * exp(-3.0 * x3)];
    self.jacobi.matrix[0][3] = [NSNumber numberWithDouble:x5 * x5];
    self.jacobi.matrix[0][4] = [NSNumber numberWithDouble:2.0 * x4 * x5];
    self.jacobi.matrix[0][5] = [NSNumber numberWithDouble:-1.0];
    self.jacobi.matrix[0][6] = [NSNumber numberWithDouble:0.0];
    self.jacobi.matrix[0][7] = [NSNumber numberWithDouble:-2.0 * cosh(2.0 * x8) * x9];
    self.jacobi.matrix[0][8] = [NSNumber numberWithDouble:-sinh(2.0 * x8)];
    self.jacobi.matrix[0][9] = [NSNumber numberWithDouble:2.0];
    self.jacobi.matrix[1][0] = [NSNumber numberWithDouble:cos(x1 * x2) * x2];
    self.jacobi.matrix[1][1] = [NSNumber numberWithDouble:cos(x1 * x2) * x1];
    self.jacobi.matrix[1][2] = [NSNumber numberWithDouble:x9 * x7];
    self.jacobi.matrix[1][3] = [NSNumber numberWithDouble:0.0];
    self.jacobi.matrix[1][4] = [NSNumber numberWithDouble:6.0 * x5];
    self.jacobi.matrix[1][5] = [NSNumber numberWithDouble:-exp(-x10 + x6) - x8 - 1.0];
    self.jacobi.matrix[1][6] = [NSNumber numberWithDouble:x3 * x9];
    self.jacobi.matrix[1][7] = [NSNumber numberWithDouble:-x6];
    self.jacobi.matrix[1][8] = [NSNumber numberWithDouble:x3 * x7];
    self.jacobi.matrix[1][9] = [NSNumber numberWithDouble:exp(-x10 + x6)];
    self.jacobi.matrix[2][0] = [NSNumber numberWithDouble:1.0];
    self.jacobi.matrix[2][1] = [NSNumber numberWithDouble:-1.0];
    self.jacobi.matrix[2][2] = [NSNumber numberWithDouble:1.0];
    self.jacobi.matrix[2][3] = [NSNumber numberWithDouble:-1.0];
    self.jacobi.matrix[2][4] = [NSNumber numberWithDouble:1.0];
    self.jacobi.matrix[2][5] = [NSNumber numberWithDouble:-1.0];
    self.jacobi.matrix[2][6] = [NSNumber numberWithDouble:1.0];
    self.jacobi.matrix[2][7] = [NSNumber numberWithDouble:-1.0];
    self.jacobi.matrix[2][8] = [NSNumber numberWithDouble:1.0];
    self.jacobi.matrix[2][9] = [NSNumber numberWithDouble:-1.0];
    self.jacobi.matrix[3][0] = [NSNumber numberWithDouble:-x5 / pow(x3 + x1, 2.0)];
    self.jacobi.matrix[3][1] = [NSNumber numberWithDouble:-2.0 * cos(x2 * x2) * x2];
    self.jacobi.matrix[3][2] = [NSNumber numberWithDouble:-x5 / pow(x3 + x1, 2.0) ];
    self.jacobi.matrix[3][3] = [NSNumber numberWithDouble:-2.0 * sin(-x9 + x4)];
    self.jacobi.matrix[3][4] = [NSNumber numberWithDouble:pow(x3 + x1, -1.0)];
    self.jacobi.matrix[3][5] = [NSNumber numberWithDouble:0.0];
    self.jacobi.matrix[3][6] = [NSNumber numberWithDouble:-2.0 * cos(x7 * x10) * sin(x7 * x10) * x10];
    self.jacobi.matrix[3][7] = [NSNumber numberWithDouble:-1.0];
    self.jacobi.matrix[3][8] = [NSNumber numberWithDouble:2.0 * sin(-x9 + x4)];
    self.jacobi.matrix[3][9] = [NSNumber numberWithDouble:-2.0 * cos(x7 * x10) * sin(x7 * x10) * x7];
    self.jacobi.matrix[4][0] = [NSNumber numberWithDouble:2.0 * x8];
    self.jacobi.matrix[4][1] = [NSNumber numberWithDouble:-2.0 * sin(x2)];
    self.jacobi.matrix[4][2] = [NSNumber numberWithDouble:2.0 * x8];
    self.jacobi.matrix[4][3] = [NSNumber numberWithDouble:pow(-x9 + x4, -2.0)];
    self.jacobi.matrix[4][4] = [NSNumber numberWithDouble:cos(x5)];
    self.jacobi.matrix[4][5] = [NSNumber numberWithDouble:x7 * exp(-x7 * (-x10 + x6))];
    self.jacobi.matrix[4][6] = [NSNumber numberWithDouble:-(x10 - x6) * exp(-x7 * (-x10 + x6))];
    self.jacobi.matrix[4][7] = [NSNumber numberWithDouble:2.0 * x3 + 2.0 * x1];
    self.jacobi.matrix[4][8] = [NSNumber numberWithDouble:-pow(-x9 + x4, -2)];
    self.jacobi.matrix[4][9] = [NSNumber numberWithDouble:-x7 * exp(-x7 * (-x10 + x6))];
    self.jacobi.matrix[5][0] = [NSNumber numberWithDouble:exp(x1 - x4 - x9)];
    self.jacobi.matrix[5][1] = [NSNumber numberWithDouble:-3.0 / 2.0 * sin(3.0 * x10 * x2) * x10];
    self.jacobi.matrix[5][2] = [NSNumber numberWithDouble:-x6];
    self.jacobi.matrix[5][3] = [NSNumber numberWithDouble:-exp(x1 - x4 - x9)];
    self.jacobi.matrix[5][4] = [NSNumber numberWithDouble:2.0 * x5 / x8];
    self.jacobi.matrix[5][5] = [NSNumber numberWithDouble:-x3];
    self.jacobi.matrix[5][6] = [NSNumber numberWithDouble:0.0];
    self.jacobi.matrix[5][7] = [NSNumber numberWithDouble:-pow(x5, 2.0) / pow(x8, 2.0)];
    self.jacobi.matrix[5][8] = [NSNumber numberWithDouble:-exp(x1 - x4 - x9)];
    self.jacobi.matrix[5][9] = [NSNumber numberWithDouble:-3.0 / 2.0 * sin(3.0 * x10 * x2) * x2];
    self.jacobi.matrix[6][0] = [NSNumber numberWithDouble:cos(x4)];
    self.jacobi.matrix[6][1] = [NSNumber numberWithDouble:3.0 * x2 * x2 * x7];
    self.jacobi.matrix[6][2] = [NSNumber numberWithDouble:1.0];
    self.jacobi.matrix[6][3] = [NSNumber numberWithDouble:-(x1 - x6) * sin(x4)];
    self.jacobi.matrix[6][4] = [NSNumber numberWithDouble:cos(x10 / x5 + x8) * x10 * pow(x5, -2.0)];
    self.jacobi.matrix[6][5] = [NSNumber numberWithDouble:-cos(x4)];
    self.jacobi.matrix[6][6] = [NSNumber numberWithDouble:pow(x2, 3.0)];
    self.jacobi.matrix[6][7] = [NSNumber numberWithDouble:-cos(x10 / x5 + x8)];
    self.jacobi.matrix[6][8] = [NSNumber numberWithDouble:0.0];
    self.jacobi.matrix[6][9] = [NSNumber numberWithDouble:-cos(x10 / x5 + x8) / x5];
    self.jacobi.matrix[7][0] = [NSNumber numberWithDouble:2.0 * x5 * (x1 - 2.0 * x6)];
    self.jacobi.matrix[7][1] = [NSNumber numberWithDouble:-x7 * exp(x2 * x7 + x10)];
    self.jacobi.matrix[7][2] = [NSNumber numberWithDouble:-2.0 * cos(-x9 + x3)];
    self.jacobi.matrix[7][3] = [NSNumber numberWithDouble:1.5];
    self.jacobi.matrix[7][4] = [NSNumber numberWithDouble:pow(x1 - 2.0 * x6, 2.0)];
    self.jacobi.matrix[7][5] = [NSNumber numberWithDouble:-4.0 * x5 * (x1 - 2.0 * x6)];
    self.jacobi.matrix[7][6] = [NSNumber numberWithDouble:-x2 * exp(x2 * x7 + x10)];
    self.jacobi.matrix[7][7] = [NSNumber numberWithDouble:0.0];
    self.jacobi.matrix[7][8] = [NSNumber numberWithDouble:2.0 * cos(-x9 + x3)];
    self.jacobi.matrix[7][9] = [NSNumber numberWithDouble:-exp(x2 * x7 + x10)];
    self.jacobi.matrix[8][0] = [NSNumber numberWithDouble:-3.0];
    self.jacobi.matrix[8][1] = [NSNumber numberWithDouble:-2.0 * x8 * x10 * x7];
    self.jacobi.matrix[8][2] = [NSNumber numberWithDouble:0.0];
    self.jacobi.matrix[8][3] = [NSNumber numberWithDouble:exp(x5 + x4)];
    self.jacobi.matrix[8][4] = [NSNumber numberWithDouble:exp(x5 + x4)];
    self.jacobi.matrix[8][5] = [NSNumber numberWithDouble:-7.0 * pow(x6, -2.0)];
    self.jacobi.matrix[8][6] = [NSNumber numberWithDouble:-2.0 * x2 * x8 * x10];
    self.jacobi.matrix[8][7] = [NSNumber numberWithDouble:-2.0 * x2 * x10 * x7];
    self.jacobi.matrix[8][8] = [NSNumber numberWithDouble:3.0];
    self.jacobi.matrix[8][9] = [NSNumber numberWithDouble:-2.0 * x2 * (x8 * x7)];
    self.jacobi.matrix[9][0] = [NSNumber numberWithDouble:x10];
    self.jacobi.matrix[9][1] = [NSNumber numberWithDouble:x9];
    self.jacobi.matrix[9][2] = [NSNumber numberWithDouble:-x8];
    self.jacobi.matrix[9][3] = [NSNumber numberWithDouble:cos(x4 + x5 + x6) * x7];
    self.jacobi.matrix[9][4] = [NSNumber numberWithDouble:cos(x4 + x5 + x6) * x7];
    self.jacobi.matrix[9][5] = [NSNumber numberWithDouble:cos(x4 + x5 + x6) * x7];
    self.jacobi.matrix[9][6] = [NSNumber numberWithDouble:sin(x4 + x5 + x6)];
    self.jacobi.matrix[9][7] = [NSNumber numberWithDouble:-x3];
    self.jacobi.matrix[9][8] = [NSNumber numberWithDouble:x2];
    self.jacobi.matrix[9][9] = [NSNumber numberWithDouble:x1];
}

// –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

- (BOOL)approachIsReached:(NSArray *)array
{
    for (NSNumber *number in array)
    {
        if (fabs([number doubleValue]) > [self.epsilon doubleValue])
        {
            return NO;
        }
    }
    return YES;
}

- (NSArray *)solveWithModified:(BOOL)isModified index:(NSUInteger)index
{
    NSDate *beginning = [NSDate date];
    
    NSMutableArray *x = [NSMutableArray arrayWithArray:self.approximation];
    NSMutableArray *prevx = [NSMutableArray arrayWithArray:self.approximation];
    NSMutableArray *solution = [NSMutableArray arrayWithArray:self.approximation];
    
    double otnos = [self max:self.approximation];
    
    BOOL recalc = NO;
    
    MKMatrix *P;
    MKMatrix *L;
    MKMatrix *U;
    
    NSUInteger k = 1;
    NSUInteger steps = 0;
    
    while (fabs(otnos) > [_epsilon doubleValue])
    {
        [self calculateFunctionInPoint:x];
        
        if (k <= index || !isModified)
        {
            [self calculateJacobiInPoint:x];
            
            printf("\n");
            
            P = [[MKMatrix alloc] init];
            L = [[MKMatrix alloc] init];
            U = [[MKMatrix alloc] init];
            
            [self.jacobi lupDecompositionWithP:P L:L U:U];
            
            recalc = NO;
        }
        
        self.previousDelta = solution;
        solution = (NSMutableArray *)[self.jacobi lupSolveWithP:P L:L U:U b:self.function];
        prevx = [x copy];
        
        for (unsigned long i = 0; i < x.count; ++i)
        {
            NSNumber *number = [NSNumber numberWithDouble:
                                ([x[i] doubleValue] + [solution[i] doubleValue])];
            [x replaceObjectAtIndex:i
                         withObject:number];
        }
        
        otnos = ([self max:x] - [self max:prevx]) / [self max:prevx];
        
        for (int i = 0; i < solution.count; ++i)
        {
            if (fabs([solution[i] doubleValue]) > fabs([self.previousDelta[i] doubleValue]))
            {
                recalc = YES;
                break;
            }
        }
        
        k++;
        steps++;
        
    }
    
    printf("Solve by \n time: %f \n steps: %ld \n \n",
           [[NSDate date] timeIntervalSinceDate:beginning], steps);
    
    return x;
}

- (double)max:(NSArray *)array
{
    double max = [array[0] doubleValue];
    
    for (NSNumber *number in array)
    {
        if (fabs([number doubleValue]) > fabs(max))
        {
            max = fabs([number doubleValue]);
        }
    }
    
    return max;
}

@end
