//
//  MKMatrix.m
//  PLU-matrix
//
//  Created by Андрей Рычков on 23.02.14.
//  Copyright (c) 2014 Andrey Rychkov. All rights reserved.
//

#import "MKMatrix.h"

@interface MKMatrix()

@property (nonatomic) NSNumber *EPS;

@end

@implementation MKMatrix

- (id)init
{
    self = [super init];
    
    if (self)
    {
        _matrix = [NSMutableArray array];
        _EPS = [NSNumber numberWithDouble:10e-9];
    }
    
    return self;
}

- (id)initWithFileName:(NSString *)fileName
{
    self = [super init];
    
    if (self)
    {
        _matrix = [[NSMutableArray alloc] initWithContentsOfFile:fileName];
        _EPS = [NSNumber numberWithDouble:10e-9];
    }
    
    return self;
}

- (void) print
{
    for (NSArray *i in _matrix)
    {
        for (NSNumber *j in i)
        {
            double d = [j doubleValue];
            printf("%.2f    ", d);
        }
        
        printf("\n");
    }
}

+ (MKMatrix *)multiplicateMatrix:(MKMatrix *)A with:(MKMatrix *)B
{
    MKMatrix *C = [[MKMatrix alloc] init];
    
    for (unsigned long i = 0; i < [A.matrix count]; ++i)
    {
        [C.matrix addObject:[NSMutableArray array]];
        for (unsigned long j = 0; j < [A.matrix count]; ++j)
        {
            double tmp = 0.0;
            for (unsigned long r = 0; r < [A.matrix count]; ++r)
            {
                tmp += [A.matrix[i][r] doubleValue] * [B.matrix[r][j] doubleValue];
            }
            [C.matrix[i] addObject:[NSNumber numberWithDouble:tmp]];
        }
    }
    
    return C;
}


// –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––



- (void)lupDecompositionWithP:(MKMatrix *)P L:(MKMatrix *)L U:(MKMatrix *)U
{
    unsigned long n = [_matrix count];
    NSMutableArray *pArray = [[NSMutableArray alloc] initWithCapacity:n];
    NSNumber *pivot;
    
    unsigned long k1 = 0;
    unsigned long i;
    
    double detP = 1.0;
    double detU = 1.0;
    
    MKMatrix *C = [[MKMatrix alloc] init];
    
    for (i = 0; i < n; ++i)
    {
        pArray[i] = [NSNumber numberWithInt:(int)i];
        
        [C.matrix addObject:[NSMutableArray arrayWithArray:[self.matrix[i] copy]]];
    }
    
    for (int k = 0; k < n; ++k)
    {
        pivot = @0.0;
        for (i = k; i < n; ++i)
        {
            if (fabs([C.matrix[i][k] doubleValue]) > pivot.doubleValue)
            {
                pivot = [NSNumber numberWithDouble:fabs([C.matrix[i][k] doubleValue])];
                k1 = i;
                
                detP = -detP;
            }
        }
        if (fabs([pivot doubleValue]) < [_EPS doubleValue])
        {
            NSLog(@"1");
            continue;
        }
        
        NSNumber *tmp = pArray[k];
        pArray[k] = pArray[k1];
        pArray[k1] = tmp;
        
        for (i = 0; i < n; ++i)
        {
            tmp = C.matrix[k][i];
            C.matrix[k][i] = C.matrix[k1][i];
            C.matrix[k1][i] = tmp;
        }
        for (i = k + 1; i < n; ++i)
        {
            C.matrix[i][k] =
            [NSNumber numberWithDouble:[C.matrix[i][k] doubleValue] /
             [C.matrix[k][k] doubleValue]];
            for (int j = k+1; j < n; ++j)
            {
                C.matrix[i][j] =
                [NSNumber numberWithDouble:[C.matrix[i][j] doubleValue] -
                 [C.matrix[i][k] doubleValue] * [C.matrix[k][j] doubleValue]];
            }
        }
    }
 
// -----------------------------------------------------------------------------
    
    for (i = 0; i < n; ++i)
    {
        [U.matrix addObject:[NSMutableArray array]];
        [L.matrix addObject:[NSMutableArray array]];
        [P.matrix addObject:[NSMutableArray array]];
        for (int j = 0; j < n; ++j)
        {
            [P.matrix[i] addObject:@0.0];
            if (i < j)
            {
                [U.matrix[i] addObject:C.matrix[i][j]];
                [L.matrix[i] addObject:@0.0];
            } else if (i > j)
            {
                [L.matrix[i] addObject:C.matrix[i][j]];
                [U.matrix[i] addObject:@0.0];
            } else if (i == j)
            {
                [U.matrix[i] addObject:C.matrix[i][j]];
                detU *= [C.matrix[i][j] doubleValue];
                [L.matrix[i] addObject:@1.0];
            }
        }
    }
    for (i = 0; i < n; ++i)
    {
        [P.matrix[i] replaceObjectAtIndex:[pArray[i] integerValue]
                                   withObject:@1.0];
    }
    
    P.detValue = [NSNumber numberWithDouble:detP];
    if ((fabs(detU)) < [_EPS doubleValue])
    {
        detU = 0.0;
    }
    U.detValue = [NSNumber numberWithDouble:detU];
    
    [P print];
    printf("\n");
    [L print];
    printf("\n");
    [U print];
    printf("\n");
}

- (void)calcDeterminantWithP:(MKMatrix *)P U:(MKMatrix *)U
{
    _detValue = [NSNumber
                     numberWithDouble:[U.detValue doubleValue] *
                     [P.detValue doubleValue]];
}

- (void)calcRankWithU:(MKMatrix *)U
{
    unsigned long rank = [U.matrix count];
    long i = 0;
    for (NSArray *array in U.matrix)
    {
        if (fabs([array[i] doubleValue]) < [_EPS doubleValue])
        {
            rank--;
        }
        ++i;
    }
    _rankValue = [NSNumber numberWithUnsignedLong:rank];
}

- (NSArray *)lupSolveWithP:(MKMatrix *)P L:(MKMatrix *)L U:(MKMatrix *)U b:(NSArray *)b
{
    unsigned long n = [_matrix count];
    NSMutableArray *x = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray *y = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray *p = [NSMutableArray arrayWithCapacity:n];
    double sum;
    
    for (unsigned long i = 0; i < n; ++i)
    {
        y[i] = [NSNumber numberWithDouble:0.0];
        x[i] = [NSNumber numberWithDouble:0.0];
        for (unsigned long j = 0; j < n; ++j)
        {
            if ([P.matrix[i][j] isEqualTo:@1.0])
            {
                p[i] = [NSNumber numberWithUnsignedLong:j];
            }
        }
    }

    if (![self.detValue isEqualTo:@0.0])
    {
        for (long i = 0; i < n; ++i)
        {
            sum = 0;
            for (long j = 0; j < i; ++j)
            {
                sum += [L.matrix[i][j] doubleValue] * [y[j] doubleValue];
            }
            sum = [b[[p[i] unsignedLongValue]] doubleValue] - sum;
            [y replaceObjectAtIndex:i
                         withObject:[NSNumber numberWithDouble:sum]];
        }
        for (long i = n - 1; i >= 0; --i)
        {
            sum = 0;
            for (long j = i + 1; j < n; ++j)
            {
                sum += [U.matrix[i][j] doubleValue] * [x[j] doubleValue];
            }
            sum = ([y[i] doubleValue] - sum) / [U.matrix[i][i] doubleValue];
            [x replaceObjectAtIndex:i
                         withObject:[NSNumber numberWithDouble:sum]];
        }
        return x;
    } else
    {
        NSMutableArray *exA = [NSMutableArray array];
        for (long i = 0; i < n; ++i)
        {
            [exA addObject:[NSMutableArray arrayWithArray:_matrix[i]]];
            [exA[i] addObject:b[i]];
        }
        
        n = [exA count];
        
        NSMutableArray *where = [NSMutableArray array];
        
        for (long i = 0; i < n; i++)
        {
            [where addObject:@-1.0];
        }
        
        for (long row = 0, col = 0; row < n && col < n; ++col, ++row)
        {
            long pivot = row;
            for (long i = row; i < n; ++i)
            {
                if (fabs([exA[i][col] doubleValue]) > fabs([exA[pivot][col] doubleValue]))
                {
                    pivot = i;
                }
            }
            if (fabs([exA[pivot][col] doubleValue]) < [_EPS doubleValue])
            {
                continue;
            }
            for (long i = col; i <= n; ++i)
            {
                NSNumber *tmp = exA[pivot][i];
                exA[pivot][i] = exA[row][i];
                exA[row][i] = tmp;
            }
            [where replaceObjectAtIndex:col
                             withObject:[NSNumber numberWithLong:row]];
            for (long i = 0; i < n; ++i)
            {
                if (i != row)
                {
                    double c = [exA[i][col] doubleValue] / [exA[row][col] doubleValue];
                    for (long j = col; j <= n; ++j)
                    {
                        [exA[i] replaceObjectAtIndex:j
                                          withObject:[NSNumber numberWithDouble:[exA[i][j] doubleValue] -
                         [exA[row][j] doubleValue] * c]];
                    }
                }
            }
        }
        
        for (long i = 0; i < n; ++i)
        {
            if (![where[i] isEqualToNumber:@-1.0])
            {
                long w = [where[i] longValue];
                double xI = [exA[w][n] doubleValue] / [exA[w][i] doubleValue];
                [x replaceObjectAtIndex:i
                             withObject:[NSNumber numberWithDouble:xI]];
            }
        }
        for (long i = 0; i < n; ++i)
        {
            sum = 0;
            for (long j = 0; j < n; ++j)
            {
                sum += [x[j] doubleValue] * [exA[i][j] doubleValue];
            }
            if (fabs(sum - [exA[i][n] doubleValue]) > [_EPS doubleValue])
            {
                NSLog(@"System is not compatible");
                return [NSArray arrayWithObjects:@0.0, nil];
            }
        }
        return x;
    }
}

- (MKMatrix *)reverseMatrixWithP:(MKMatrix *)P L:(MKMatrix *)L U:(MKMatrix *)U
{
    if ([self.detValue isEqualToNumber:@0.0])
    {
        NSLog(@"Matrix is singular");
        return [[MKMatrix alloc] init];
    }
    
    MKMatrix *I = [[MKMatrix alloc] init];
    long n = [_matrix count];
    NSMutableArray *b = [NSMutableArray array];
    
    [b addObject:@1.0];
    [I.matrix addObject:[NSMutableArray array]];
    for (long i = 1; i < n; ++i)
    {
        [b addObject:@0.0];
        [I.matrix addObject:[NSMutableArray array]];
    }
    
    for (long i = 0; i < n; ++i)
    {
        NSArray *col = [self lupSolveWithP:P L:L U:U b:b];
        
        for (long j = 0; j < n; ++j)
        {
            [I.matrix[j] addObject:col[j]];
        }
        
        if (i + 1 < n)
        {
            NSNumber *tmp = b[i];
            b[i] = b[i+1];
            b[i+1] = tmp;
        }
    }
    
    return I;
}

- (NSNumber *)calcConditionNumber:(MKMatrix *)I
{
    double condA = 0;
    double condI = 0;
    long n = [_matrix count];
    
    for (long i = 0; i < n; ++i)
    {
        for (long j = 0; j < n; ++j)
        {
            condA += [_matrix[i][j] doubleValue] * [_matrix[i][j] doubleValue];
            condI += [I.matrix[i][j] doubleValue] * [I.matrix[i][j] doubleValue];
        }
    }
    
    condA = sqrt(condA);
    condI = sqrt(condI);
    
    return [NSNumber numberWithDouble:condA * condI];
}

@end
