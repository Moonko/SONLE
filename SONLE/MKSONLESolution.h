//
//  MKSONLESolution.h
//  SONLE
//
//  Created by Андрей Рычков on 16.03.14.
//  Copyright (c) 2014 Andrey Rychkov. All rights reserved.
//

#import <Foundation/Foundation.h>

@class MKMatrix;

@interface MKSONLESolution : NSObject

@property (nonatomic) NSNumber *epsilon;
@property (nonatomic) NSArray *approximation;
@property (nonatomic) NSMutableArray *function;
@property (nonatomic) NSMutableArray *previousDelta;
@property (nonatomic) MKMatrix *jacobi;

- (id)initWithEpsilon:(NSNumber *)eps
        approximation:(NSArray *)approximation
    numberOfVariables:(unsigned long)number;

- (NSArray *)solveWithModified:(BOOL)isModified index:(NSUInteger)index;

@end
