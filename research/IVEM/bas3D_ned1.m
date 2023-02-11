function bas_curl = bas3D_ned1(bas)

nt = size(bas,1);  
bas_curl_x = zeros(nt,4,6); bas_curl_y = zeros(nt,4,6); bas_curl_z = zeros(nt,4,6);
e_ind = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];

for j = 1:6
    for i = 1:4
        bas_curl_x(:,i,j) = bas(:,i,e_ind(j,1)).*bas(:,2,e_ind(j,2)) -...
            bas(:,i,e_ind(j,2)).*bas(:,2,e_ind(j,1));
        bas_curl_y(:,i,j) = bas(:,i,e_ind(j,1)).*bas(:,3,e_ind(j,2)) -...
            bas(:,i,e_ind(j,2)).*bas(:,3,e_ind(j,1));
        bas_curl_z(:,i,j) = bas(:,i,e_ind(j,1)).*bas(:,4,e_ind(j,2)) -...
            bas(:,i,e_ind(j,2)).*bas(:,4,e_ind(j,1));
    end
end

bas_curl.bas_curl_x=bas_curl_x;
bas_curl.bas_curl_y=bas_curl_y;
bas_curl.bas_curl_z=bas_curl_z;

