//
//  MKMatrix.h
//  PLU-matrix
//
//  Created by Андрей Рычков on 23.02.14.
//  Copyright (c) 2014 Andrey Rychkov. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface MKMatrix : NSObject

@property (nonatomic) NSMutableArray *matrix;
@property (nonatomic) NSNumber *detValue;
@property (nonatomic) NSNumber *rankValue;

- (id)initWithFileName:(NSString *)fileName;

- (void)print;

+ (MKMatrix *)multiplicateMatrix:(MKMatrix *)A with:(MKMatrix *)B;

- (void)lupDecompositionWithP:(MKMatrix *)P L:(MKMatrix *)L U:(MKMatrix *)U;

- (void)calcDeterminantWithP:(MKMatrix *)P U:(MKMatrix *)U;

- (void)calcRankWithU:(MKMatrix *)U;

- (NSArray *)lupSolveWithP:(MKMatrix *)P L:(MKMatrix *)L U:(MKMatrix *)U b:(NSArray *)b;

- (MKMatrix *)reverseMatrixWithP:(MKMatrix *)P L:(MKMatrix *)L U:(MKMatrix *)U;

- (NSNumber *)calcConditionNumber:(MKMatrix *)I;

@end
