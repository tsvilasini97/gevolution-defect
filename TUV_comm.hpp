#ifndef TUV_COMM_HPP
#define TUV_COMM_HPP

void TuvProjectionCICNGP_comm(Field<Real> * vel)
{

    if(vel->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSend;
    Real *bufferRec;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = vel->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=vel->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int iref;
    int imax;




    int comp=10;


    iref =sizeLocalGross[0] - halo ;
    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int c=0;c<comp;c++)(*vel)(setIndex(sizeLocalGross,halo,j,k),c) += (*vel)(setIndex(sizeLocalGross,iref,j,k),c);
        }
    }


    //send halo in direction Y
    bufferSizeY =  (long)(sizeLocalOne[2]-1)*sizeLocal[0] * comp;
    bufferSizeZ = sizeLocal[0] * sizeLocal[1] * comp;
    if(bufferSizeY>bufferSizeZ)
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY);
    }
    else
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ);
    }


    //pack data
    imax=sizeLocalGross[0]-2* halo;
    iref=sizeLocalGross[1]- halo;
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+k*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
        }
    }


    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);


    //unpack data
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*vel)(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[c+comp*(i+k*imax)];
        }

    }

    //send halo in direction Z

    //pack data

    //cout<<"okok"<<endl;

    iref=sizeLocalGross[2]-halo;
    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+j*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
        }
    }

    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);


    //unpack data

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[c+comp*(i+j*imax)];
        }
    }

    free(bufferRec);
    free(bufferSend);


}

void TUV_comm(Field<Real> * Tij)
{

    if(Tij->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSend;
    Real *bufferRec;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = Tij->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=Tij->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int iref;
    int imax;




    int comp=10;

    iref = sizeLocalGross[0]-halo;
    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,halo,j,k),c) += (*Tij)(setIndex(sizeLocalGross,iref,j,k),c);
        }
    }


    //send halo in direction Y
    bufferSizeY =  (long)(sizeLocalOne[2]-1)*sizeLocal[0] * comp;
    bufferSizeZ = sizeLocal[0] * sizeLocal[1] * comp;
    if(bufferSizeY>bufferSizeZ)
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY);
    }
    else
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ);
    }


    //pack data
    imax=sizeLocalGross[0]-2* halo;
    iref=sizeLocalGross[1]- halo;
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
        }
    }


    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);


    //unpack data
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[c+comp*(i+k*imax)];
        }

    }

//    //send halo in direction Z

    //pack data
    iref=sizeLocalGross[2]-halo;
    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
        }
    }

    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);


    //unpack data

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[c+comp*(i+j*imax)];
        }
    }

    free(bufferRec);
    free(bufferSend);



}


#endif
